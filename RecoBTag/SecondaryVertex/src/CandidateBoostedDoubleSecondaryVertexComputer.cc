#include "RecoBTag/SecondaryVertex/interface/CandidateBoostedDoubleSecondaryVertexComputer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/BTauGenericMVAJetTagComputerRcd.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


CandidateBoostedDoubleSecondaryVertexComputer::CandidateBoostedDoubleSecondaryVertexComputer(const edm::ParameterSet & parameters) :
  beta_(parameters.getParameter<double>("beta")),
  R0_(parameters.getParameter<double>("R0")),
  njettiness_(fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta_,R0_)),
  maxSVDeltaRToJet_(parameters.getParameter<double>("maxSVDeltaRToJet")),
  useCondDB_(parameters.getParameter<bool>("useCondDB")),
  gbrForestLabel_(parameters.existsAs<std::string>("gbrForestLabel") ? parameters.getParameter<std::string>("gbrForestLabel") : ""),
  weightFile_(parameters.existsAs<edm::FileInPath>("weightFile") ? parameters.getParameter<edm::FileInPath>("weightFile") : edm::FileInPath()),
  useGBRForest_(parameters.existsAs<bool>("useGBRForest") ? parameters.getParameter<bool>("useGBRForest") : false),
  useAdaBoost_(parameters.existsAs<bool>("useAdaBoost") ? parameters.getParameter<bool>("useAdaBoost") : false)
{
  uses(0, "ipTagInfos");
  uses(1, "svTagInfos");

  mvaID.reset(new TMVAEvaluator());
}

void CandidateBoostedDoubleSecondaryVertexComputer::initialize(const JetTagComputerRecord & record)
{
  // variable names and order need to be the same as in the training
  std::vector<std::string> variables({"z_ratio",
                                      "trackSipdSig_3","trackSipdSig_2","trackSipdSig_1","trackSipdSig_0",
                                      "trackSipdSig_1_0","trackSipdSig_0_0","trackSipdSig_1_1","trackSipdSig_0_1",
                                      "trackSip2dSigAboveCharm_0","trackSip2dSigAboveBottom_0","trackSip2dSigAboveBottom_1",
                                      "tau0_trackEtaRel_0","tau0_trackEtaRel_1","tau0_trackEtaRel_2",
                                      "tau1_trackEtaRel_0","tau1_trackEtaRel_1","tau1_trackEtaRel_2",
                                      "tau_vertexMass_0","tau_vertexEnergyRatio_0","tau_vertexDeltaR_0","tau_flightDistance2dSig_0",
                                      "tau_vertexMass_1","tau_vertexEnergyRatio_1","tau_flightDistance2dSig_1",
                                      "jetNTracks","nSV"});
  // book TMVA readers
  std::vector<std::string> spectators({"massPruned", "flavour", "nbHadrons", "ptPruned", "etaPruned"});

  if (useCondDB_)
  {
     const GBRWrapperRcd & gbrWrapperRecord = record.getRecord<GBRWrapperRcd>();

     edm::ESHandle<GBRForest> gbrForestHandle;
     gbrWrapperRecord.get(gbrForestLabel_.c_str(), gbrForestHandle);

     mvaID->initializeGBRForest(gbrForestHandle.product(), variables, spectators, useAdaBoost_);
  }
  else
    mvaID->initialize("Color:Silent:Error", "BDT", weightFile_.fullPath(), variables, spectators, useGBRForest_, useAdaBoost_);
}

float CandidateBoostedDoubleSecondaryVertexComputer::discriminator(const TagInfoHelper & tagInfo) const
{
  // get TagInfos
//   const reco::CandIPTagInfo              & ipTagInfo = tagInfo.get<reco::CandIPTagInfo>(0);
//   const reco::CandSecondaryVertexTagInfo & svTagInfo = tagInfo.get<reco::CandSecondaryVertexTagInfo>(1);

  // default discriminator value
  float value = -10.;

  // default variable values
//   float z_ratio = -3.;
//   float trackSip3dSig_3 = -50., trackSip3dSig_2 = -50., trackSip3dSig_1 = -50., trackSip3dSig_0 = -50.;
//   float tau2_trackSip3dSig_0 = -50., tau1_trackSip3dSig_0 = -50., tau2_trackSip3dSig_1 = -50., tau1_trackSip3dSig_1 = -50.;
//   float trackSip2dSigAboveCharm_0 = -19., trackSip2dSigAboveBottom_0 = -19., trackSip2dSigAboveBottom_1 = -19.;
//   float tau1_trackEtaRel_0 = -1., tau1_trackEtaRel_1 = -1., tau1_trackEtaRel_2 = -1.;
//   float tau2_trackEtaRel_0 = -1., tau2_trackEtaRel_1 = -1., tau2_trackEtaRel_2 = -1.;
//   float tau1_vertexMass = -1., tau1_vertexEnergyRatio = -1., tau1_vertexDeltaR = -1., tau1_flightDistance2dSig = -1.;
//   float tau_vertexMass_1 = -1., tau_vertexEnergyRatio_1 = -1., tau_flightDistance2dSig_1 = -1.;
//   float jetNTracks = 0, nSV = 0;

  // get the jet reference
//   const reco::JetBaseRef jet = svTagInfo.jet();

//   std::vector<fastjet::PseudoJet> currentAxes;
//   float tau2, tau1;
//   // calculate N-subjettiness
//   calcNsubjettiness(jet, tau1, tau2, currentAxes);
//   if (tau1 != 0.) tau21 = tau2/tau1;
// 
//   const std::vector<reco::CandidatePtr> & selectedTracks( ipTagInfo.selectedTracks() );
//   size_t trackSize = selectedTracks.size();
//   const reco::VertexRef & vertexRef = ipTagInfo.primaryVertex();
//   reco::TrackKinematics allKinematics;
// 
//   for (size_t itt=0; itt < trackSize; ++itt)
//   {
//     const reco::Track & ptrack = *(reco::btag::toTrack(selectedTracks[itt]));
//     const reco::CandidatePtr ptrackRef = selectedTracks[itt];
// 
//     float track_PVweight = 0.;
//     setTracksPV(ptrackRef, vertexRef, track_PVweight);
//     if (track_PVweight>0.) { allKinematics.add(ptrack, track_PVweight); }
//   }
// 
//   math::XYZVector jetDir = jet->momentum().Unit();
// 
//   std::map<double, size_t> VTXmass;
//   for (size_t vtx = 0; vtx < svTagInfo.nVertices(); ++vtx)
//   {
//     vertexNTracks += (svTagInfo.secondaryVertex(vtx)).numberOfSourceCandidatePtrs();
//     GlobalVector flightDir = svTagInfo.flightDirection(vtx);
//     if (reco::deltaR2(flightDir, jetDir)<(maxSVDeltaRToJet_*maxSVDeltaRToJet_))
//     {
//       ++contSV;
//       VTXmass[svTagInfo.secondaryVertex(vtx).p4().mass()]=vtx;
//     }
//   }
// 
//   int cont=0;
//   GlobalVector flightDir_0, flightDir_1;
//   reco::Candidate::LorentzVector SV_p4_0 , SV_p4_1;
//   for ( std::map<double, size_t>::reverse_iterator iVtx=VTXmass.rbegin(); iVtx!=VTXmass.rend(); ++iVtx)
//   {
//     ++cont;
//     const reco::VertexCompositePtrCandidate &vertex = svTagInfo.secondaryVertex(iVtx->second);
//     reco::TrackKinematics vtxKinematics;
//     vertexKinematics(vertex, vtxKinematics);
//     math::XYZTLorentzVector allSum = allKinematics.weightedVectorSum();
//     math::XYZTLorentzVector vertexSum = vtxKinematics.weightedVectorSum();
//     if (cont==1)
//     {
//       SV_mass_0 = vertex.p4().mass()  ;
//       SV_EnergyRatio_0 = vertexSum.E() / allSum.E();
//       SV_pt_0 = vertex.p4().pt();
//       flightDir_0 = svTagInfo.flightDirection(iVtx->second);
//       SV_p4_0 = vertex.p4();
// 
//       if (reco::deltaR2(flightDir_0,currentAxes[1])<reco::deltaR2(flightDir_0,currentAxes[0]))
//         tau_dot = (currentAxes[1].px()*flightDir_0.x()+currentAxes[1].py()*flightDir_0.y()+currentAxes[1].pz()*flightDir_0.z())/(sqrt(currentAxes[1].modp2())*flightDir_0.mag());
//       else
//         tau_dot = (currentAxes[0].px()*flightDir_0.x()+currentAxes[0].py()*flightDir_0.y()+currentAxes[0].pz()*flightDir_0.z())/(sqrt(currentAxes[0].modp2())*flightDir_0.mag());
//     }
//     if (cont==2)
//     {
//       SV_EnergyRatio_1= vertexSum.E() / allSum.E();
//       flightDir_1 = svTagInfo.flightDirection(iVtx->second);
//       SV_p4_1 = vertex.p4();
//       z_ratio = reco::deltaR(flightDir_0,flightDir_1)*SV_pt_0/(SV_p4_0+SV_p4_1).mass();
//       break;
//     }
//   }

//   std::map<std::string,float> inputs;
//   inputs["z_ratio"] = z_ratio;
//   inputs["trackSipdSig_3"] = trackSip3dSig_3;
//   inputs["trackSipdSig_2"] = trackSip3dSig_2;
//   inputs["trackSipdSig_1"] = trackSip3dSig_1;
//   inputs["trackSipdSig_0"] = trackSip3dSig_0;
//   inputs["trackSipdSig_1_0"] = tau2_trackSip3dSig_0;
//   inputs["trackSipdSig_0_0"] = tau1_trackSip3dSig_0;
//   inputs["trackSipdSig_1_1"] = tau2_trackSip3dSig_1;
//   inputs["trackSipdSig_0_1"] = tau1_trackSip3dSig_1;
//   inputs["trackSip2dSigAboveCharm_0"] = trackSip2dSigAboveCharm_0;
//   inputs["trackSip2dSigAboveBottom_0"] = trackSip2dSigAboveBottom_0;
//   inputs["trackSip2dSigAboveBottom_1"] = trackSip2dSigAboveBottom_1;
//   inputs["tau1_trackEtaRel_0"] = tau2_trackEtaRel_0;
//   inputs["tau1_trackEtaRel_1"] = tau2_trackEtaRel_1;
//   inputs["tau1_trackEtaRel_2"] = tau2_trackEtaRel_2;
//   inputs["tau0_trackEtaRel_0"] = tau1_trackEtaRel_0;
//   inputs["tau0_trackEtaRel_1"] = tau1_trackEtaRel_1;
//   inputs["tau0_trackEtaRel_2"] = tau1_trackEtaRel_2;
//   inputs["tau_vertexMass_0"] = tau1_vertexMass;
//   inputs["tau_vertexEnergyRatio_0"] = tau1_vertexEnergyRatio;
//   inputs["tau_vertexDeltaR_0"] = tau1_vertexDeltaR;
//   inputs["tau_flightDistance2dSig_0"] = tau1_flightDistance2dSig;
//   inputs["tau_vertexMass_1"] = tau2_vertexMass;
//   inputs["tau_vertexEnergyRatio_1"] = tau2_vertexEnergyRatio;
//   inputs["tau_flightDistance2dSig_1"] = tau2_flightDistance2dSig;
//   inputs["jetNTracks"] = jetNTracks;
//   inputs["nSV"] = nSV;

  // evaluate the MVA
//   value = mvaID->evaluate(inputs);

  // return the final discriminator value
  return value;
}


void CandidateBoostedDoubleSecondaryVertexComputer::calcNsubjettiness(const reco::JetBaseRef & jet, float & tau1, float & tau2, std::vector<fastjet::PseudoJet> & currentAxes) const
{
  std::vector<fastjet::PseudoJet> fjParticles;

  // loop over jet constituents and push them in the vector of FastJet constituents
  for(const reco::CandidatePtr & daughter : jet->daughterPtrVector())
  {
    if ( daughter.isNonnull() && daughter.isAvailable() )
      fjParticles.push_back( fastjet::PseudoJet( daughter->px(), daughter->py(), daughter->pz(), daughter->energy() ) );
    else
      edm::LogWarning("MissingJetConstituent") << "Jet constituent required for N-subjettiness computation is missing!";
  }

  // calculate N-subjettiness
  tau1 = njettiness_.getTau(1, fjParticles);
  tau2 = njettiness_.getTau(2, fjParticles);
  currentAxes = njettiness_.currentAxes();
}


void CandidateBoostedDoubleSecondaryVertexComputer::setTracksPVBase(const reco::TrackRef & trackRef, const reco::VertexRef & vertexRef, float & PVweight) const
{
  PVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::Vertex::trackRef_iterator IT;

  const reco::Vertex & vtx = *(vertexRef);
  // loop over tracks in vertices
  for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
  {
    const reco::TrackBaseRef & baseRef = *it;
    // one of the tracks in the vertex is the same as the track considered in the function
    if( baseRef == trackBaseRef )
    {
      PVweight = vtx.trackWeight(baseRef);
      break;
    }
  }
}


void CandidateBoostedDoubleSecondaryVertexComputer::setTracksPV(const reco::CandidatePtr & trackRef, const reco::VertexRef & vertexRef, float & PVweight) const
{
  PVweight = 0.;

  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(trackRef.get());

  if(pcand) // MiniAOD case
  {
    if( pcand->fromPV() == pat::PackedCandidate::PVUsedInFit )
    {
      PVweight = 1.;
    }
  }
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(trackRef.get());

    setTracksPVBase(pfcand->trackRef(), vertexRef, PVweight);
  }
}


void CandidateBoostedDoubleSecondaryVertexComputer::vertexKinematics(const reco::VertexCompositePtrCandidate & vertex, reco::TrackKinematics & vtxKinematics) const
{
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
    const reco::Track& mytrack = *(*track)->bestTrack();
    vtxKinematics.add(mytrack, 1.0);
  }
}
