#ifndef MUON_PAIR_MAKER_H
#define MUON_PAIR_MAKER_H

#include "TreeAnalyzer.h"
#include "HistoBins.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"
#include "FemtoDstFormat/FemtoBTofPidTraits.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

// Analyzers
// #include "MuonPairMaker/PairHistogramMaker.h"
// #include "MuonPairMaker/TrackHistogramMaker.h"
// #include "MuonPairMaker/MtdHistogramMaker.h"

#include "Filters/TrackFilter.h"
#include "Filters/LowPtMuonFilter.h"
#include "Filters/MuonBDTFilter.h"
#include "Filters/MuonMLPFilter.h"
#include "Filters/MtdFilter.h"

#include <map>

#include "TLorentzVector.h"
#include "TNTuple.h"

#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"

class MuonPairMaker : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoBTofPidTraits> _rBTofPid;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	TrackFilter _lowPtTrackFilter;
	TrackFilter _mtdTrackFilter;
	LowPtMuonFilter _lowFilter;
	MuonBDTFilter _bdtFilter;
	MuonMLPFilter _mlpFilter;
	MtdFilter _mtdFilter;

	map<int, int> cMap;

	TNtuple * tNtuple = nullptr;

	bool makeHistograms = false;

public:
	virtual const char* classname() const {return "MuonPairMaker";}
	MuonPairMaker() {}
	~MuonPairMaker() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rHelices.setup( chain, "Helices" );
		_rBTofPid.setup( chain, "BTofPidTraits" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		cMap = config.getIntMap( nodePath + ".CentralityMap" );
		for ( int i = 0; i < 16; i++ )
			LOG_F( INFO, "cMap[%d] = %d", i, cMap[i] );

		_lowPtTrackFilter.load( config, nodePath + ".LowPtTrackFilter" );
		_mtdTrackFilter.load( config, nodePath + ".MtdTrackFilter" );
		_lowFilter.load( config, nodePath + ".LowPtMuonFilter" );
		_bdtFilter.load( config, nodePath + ".MuonBDTFilter" );
		_mlpFilter.load( config, nodePath + ".MuonMLPFilter" );
		_mtdFilter.load( config, nodePath + ".MtdPidFilter" );
		book->cd();

		if ( config.getBool( nodePath + ".output:tuple", false ) )
			tNtuple = new TNtuple( "Pairs", "", "lh1:lh2:tof1:tof2:mass:pT:cs" );
		makeHistograms = config.getBool( nodePath + ".output:histograms", true );
	}


protected:

	vector<FemtoTrackProxy> pos_mtd;
	vector<FemtoTrackProxy> neg_mtd;
	vector<FemtoTrackProxy> pos_smtd;
	vector<FemtoTrackProxy> neg_smtd;
	vector<FemtoTrackProxy> pos_tof;
	vector<FemtoTrackProxy> neg_tof;

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();

	}

	virtual void makePairs( vector<FemtoTrackProxy> &col, string prefix ){
		TLorentzVector lv1, lv2, lv;

		for ( size_t i = 0; i < col.size(); i++ ){
			for ( size_t j = i; j < col.size(); j++ ){
				if ( i == j ) continue;
				FemtoTrackProxy _p1 = col[i];
				FemtoTrackProxy _p2 = col[j];

				lv1.SetPtEtaPhiM( _p1._track->mPt, _p1._track->mEta, _p1._track->mPhi, 0.105 );
				lv2.SetPtEtaPhiM( _p2._track->mPt, _p2._track->mEta, _p2._track->mPhi, 0.105 );
				lv = lv1 + lv2;

				if (makeHistograms) {
					book->fill( prefix + "_pt_mass", lv.M(), lv.Pt() );
					book->fill( prefix + "_lh_mass", lv.M(), _p1._lh, _p2._lh );
				}
				
				if ( nullptr != tNtuple ){
					float tof1 = -999;
					float tof2 = -999;
					if ( nullptr!= _p1._mtdPid )
						tof1 = _p1._mtdPid->mDeltaTimeOfFlight;
					if ( nullptr!= _p2._mtdPid )
						tof2 = _p2._mtdPid->mDeltaTimeOfFlight;
					float cs = _p1._track->charge() + _p2._track->charge();
					
					tNtuple->Fill( _p1._lh, _p2._lh, tof1, tof2, lv.M(), lv.Pt(), cs );
				}
				
			} // j	
		}// i
	}


	virtual void makePairs( vector<FemtoTrackProxy> &col1, vector<FemtoTrackProxy> &col2, string prefix ){
		TLorentzVector lv1, lv2, lv;
		for ( FemtoTrackProxy& _proxy1 : col1 ){
			for ( FemtoTrackProxy& _proxy2 : col2 ){
				if ( fabs( _proxy1._track->mPt - _proxy2._track->mPt ) < 0.001 && fabs( _proxy1._track->mEta - _proxy2._track->mEta ) < 0.001 && fabs( _proxy1._track->mPhi - _proxy2._track->mPhi ) < 0.001 )
					continue;
				lv1.SetPtEtaPhiM( _proxy1._track->mPt, _proxy1._track->mEta, _proxy1._track->mPhi, 0.105 );
				lv2.SetPtEtaPhiM( _proxy2._track->mPt, _proxy2._track->mEta, _proxy2._track->mPhi, 0.105 );
				lv = lv1 + lv2;


				if (makeHistograms) {
					book->fill( prefix + "_pt_mass", lv.M(), lv.Pt() );
					book->fill( prefix + "_lh_mass", lv.M(), _proxy1._lh, _proxy2._lh );
				}
				if ( nullptr != tNtuple ){
					float tof1 = -999;
					float tof2 = -999;
					if ( nullptr!= _proxy1._mtdPid )
						tof1 = _proxy1._mtdPid->mDeltaTimeOfFlight;
					if ( nullptr!= _proxy2._mtdPid )
						tof2 = _proxy2._mtdPid->mDeltaTimeOfFlight;
					float cs = _proxy1._track->charge() + _proxy2._track->charge();
					tNtuple->Fill( _proxy1._lh, _proxy2._lh, tof1, tof2, lv.M(), lv.Pt(), cs );
				}
			} // loop col1
		} // loop col2
	}

	virtual void analyzeEvent(){
	
		_event = _rEvent.get();
		book->fill( "Events", 1 );

		
		if ( cMap.count( _event->mBin16 ) == 0 ) return;
		
		int mappedCen = cMap[ _event->mBin16 ];
		
		// LOG_F( INFO, "cMap[ %d ] = %d", _event->mBin16, mappedCen );
		book->fill( "mBin16", _event->mBin16 );
		book->fill( "mMappedCen", mappedCen );
		
		if ( mappedCen < 0 ) return;

		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;

		size_t nTOF = 0;
		size_t nMTD = 0;



		pos_smtd.clear();
		neg_smtd.clear();
		pos_mtd.clear();
		neg_mtd.clear();
		pos_tof.clear();
		neg_tof.clear();

		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rBTofPid );
			_proxy.setMtdPidTraits( _rMtdPid );
			_proxy._lh = -1;

			int charge = _proxy._track->charge();
			
			if ( _proxy._track->mBTofPidTraitsIndex >= 0 && _lowPtTrackFilter.pass( _proxy ) ){
				double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
				double zb = _lowFilter.zb( _proxy, "mu" );
				
				if (makeHistograms) book->fill( "zb_p", p, zb );
				if ( _lowFilter.pass( _proxy ) ){
					if (makeHistograms) book->fill( "zb_p_signal", p,zb );
					if ( charge > 0 )
						pos_tof.push_back( _proxy );
					else 
						neg_tof.push_back( _proxy );
					nTOF++;
				}
			}
			
			if ( _mtdTrackFilter.pass( _proxy ) &&  _proxy._track->mMtdPidTraitsIndex >= 0 && _proxy._mtdPid->mTriggerFlag >= 0 ){
				// DeltaZ > 60 is HACK until new MC is trained
				float bdt = _bdtFilter.evaluate( _proxy );
				float mlp = _mlpFilter.evaluate( _proxy );
				if (makeHistograms) {
					book->fill( "BDT", bdt );
					book->fill( "MLP", mlp );
					book->fill( "MLP_vs_BDT", bdt, mlp );

					book->fill( "dY_BDT", bdt, _proxy._mtdPid->mDeltaY );
					book->fill( "dZ_BDT", bdt, _proxy._mtdPid->mDeltaZ );
					book->fill( "dTOF_BDT", bdt, _proxy._mtdPid->mDeltaTimeOfFlight );
					book->fill( "nSigmaPi_BDT", bdt, _proxy._track->nSigmaPion() );
					book->fill( "nHitsFit_BDT", bdt, fabs(_proxy._track->mNHitsFit) );
					book->fill( "dca_BDT", bdt, _proxy._track->gDCA() );
					book->fill( "cell_BDT", bdt, _proxy._mtdPid->cell() );
					book->fill( "module_BDT", bdt, _proxy._mtdPid->module() );
					book->fill( "backleg_BDT", bdt, _proxy._mtdPid->backleg() );
					book->fill( "pT_BDT", bdt, _proxy._track->mPt );
					book->fill( "charge_BDT", bdt, _proxy._track->charge() );

					book->fill( "dY_MLP", mlp, _proxy._mtdPid->mDeltaY );
					book->fill( "dZ_MLP", mlp, _proxy._mtdPid->mDeltaZ );
					book->fill( "dTOF_MLP", mlp, _proxy._mtdPid->mDeltaTimeOfFlight );
					book->fill( "nSigmaPi_MLP", mlp, _proxy._track->nSigmaPion() );
					book->fill( "nHitsFit_MLP", mlp, fabs(_proxy._track->mNHitsFit) );
					book->fill( "dca_MLP", mlp, _proxy._track->gDCA() );
					book->fill( "cell_MLP", mlp, _proxy._mtdPid->cell() );
					book->fill( "module_MLP", mlp, _proxy._mtdPid->module() );
					book->fill( "backleg_MLP", mlp, _proxy._mtdPid->backleg() );
					book->fill( "pT_MLP", mlp, _proxy._track->mPt );
					book->fill( "charge_MLP", mlp, _proxy._track->charge() );
				}
				_proxy._lh = mlp;
				if ( _mlpFilter.pass( _proxy ) ){
					if ( charge > 0 )
						pos_mtd.push_back( _proxy );
					else 
						neg_mtd.push_back( _proxy );
					nMTD++;
				}
				
			}

			if ( _mtdTrackFilter.pass( _proxy ) && _mtdFilter.pass( _proxy ) ){
				if ( charge > 0 )
					pos_smtd.push_back( _proxy );
				else 
					neg_smtd.push_back( _proxy );
			}


		} // loop on tracks

		book->fill( "nTof_vs_nMtd", nMTD, nTOF );
		book->fill( "tof_pos_vs_neg", neg_tof.size(), pos_tof.size() );
		book->fill( "mtd_pos_vs_neg", neg_mtd.size(), pos_mtd.size() );

		makePairs( pos_tof, neg_mtd, "uls" );
		makePairs( neg_tof, pos_mtd, "uls" );
		makePairs( pos_tof, pos_mtd, "ls" );
		makePairs( neg_tof, neg_mtd, "ls" );

		// makePairs( pos_tof, neg_tof, "tof_uls" );
		// makePairs( neg_tof, "tof_ls" );
		// makePairs( pos_tof, "tof_ls" );

		makePairs( pos_mtd, neg_mtd, "mtd_uls" );
		makePairs( neg_mtd, "mtd_ls" );
		makePairs( pos_mtd, "mtd_ls" );

		makePairs( pos_smtd, neg_smtd, "smtd_uls" );
		makePairs( neg_smtd, "smtd_ls" );
		makePairs( pos_smtd, "smtd_ls" );


	}


	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}

		book->cd();
		tNtuple->Write();
	}
	
};

#endif
