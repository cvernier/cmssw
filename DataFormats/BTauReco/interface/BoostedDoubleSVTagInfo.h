#include <string>
#include <map>

#ifndef RecoBTag_SecondaryVertex_interface_CandidateBoostedDoubleSecondaryVertexComputer_h
#define RecoBTag_SecondaryVertex_interface_CandidateBoostedDoubleSecondaryVertexComputer_h

#define DECLARE_VARIABLE(NAME,TYPE)		\
	private: \
	TYPE NAME ## _;			\
	public: \
	const TYPE & NAME() const { return NAME ## _; } \
	void NAME(const TYPE val) { NAME ## _ = val; } 

// ----------------------------------------------------------------------------------------------------
class StoredBoostedDoubleSVTagInfo {
public:
	friend class CandidateBoostedDoubleSecondaryVertexComputer;

	StoredBoostedDoubleSVTagInfo();
	~StoredBoostedDoubleSVTagInfo(); 

	DECLARE_VARIABLE(mva      ,float);

	DECLARE_VARIABLE(z_ratio, float);
	DECLARE_VARIABLE(trackSipdSig_3, float);
	DECLARE_VARIABLE(trackSipdSig_2, float);
	DECLARE_VARIABLE(trackSipdSig_1, float);
	DECLARE_VARIABLE(trackSipdSig_0, float) ;
	DECLARE_VARIABLE(trackSipdSig_1_0, float) ;
	DECLARE_VARIABLE(trackSipdSig_0_0, float) ;
	DECLARE_VARIABLE(trackSipdSig_1_1, float) ;
	DECLARE_VARIABLE(trackSipdSig_0_1, float) ;
	DECLARE_VARIABLE(trackSip2dSigAboveCharm_0, float);
	DECLARE_VARIABLE(trackSip2dSigAboveBottom_0, float);
	DECLARE_VARIABLE(trackSip2dSigAboveBottom_1, float);
	DECLARE_VARIABLE(tau1_trackEtaRel_0, float);
	DECLARE_VARIABLE(tau1_trackEtaRel_1, float);
	DECLARE_VARIABLE(tau1_trackEtaRel_2, float) ;
	DECLARE_VARIABLE(tau0_trackEtaRel_0, float) ;
	DECLARE_VARIABLE(tau0_trackEtaRel_1, float) ;
	DECLARE_VARIABLE(tau0_trackEtaRel_2, float) ;
	DECLARE_VARIABLE(tau_vertexMass_0, float) ;
	DECLARE_VARIABLE(tau_vertexEnergyRatio_0, float) ;
	DECLARE_VARIABLE(tau_vertexDeltaR_0, float) ;
	DECLARE_VARIABLE(tau_flightDistance2dSig_0, float) ;
	DECLARE_VARIABLE(tau_vertexMass_1, float) ;
	DECLARE_VARIABLE(tau_vertexEnergyRatio_1, float) ; 
	DECLARE_VARIABLE(tau_flightDistance2dSig_1, float) ;
	DECLARE_VARIABLE(jetNTracks, int) ;
	DECLARE_VARIABLE(nSV, int) ;

        DECLARE_VARIABLE(massPruned,float);
        DECLARE_VARIABLE(flavour, float);
        DECLARE_VARIABLE(nbHadrons, int);
        DECLARE_VARIABLE(ptPruned, float);
        DECLARE_VARIABLE(etaPruned, float);

};

// ----------------------------------------------------------------------------------------------------
class BoostedDoubleSVTagInfo : public StoredBoostedDoubleSVTagInfo {
	public:
		friend class CandidateBoostedDoubleSecondaryVertexComputer;

		BoostedDoubleSVTagInfo();
		~BoostedDoubleSVTagInfo(); 

		BoostedDoubleSVTagInfo & operator= (const StoredBoostedDoubleSVTagInfo & lhs) { ((StoredBoostedDoubleSVTagInfo &)(*this)) = lhs; return *this;}


//		DECLARE_VARIABLE(jetPhi   ,float);
//		DECLARE_VARIABLE(jetM     ,float);
};

#undef DECLARE_VARIABLE

#endif
