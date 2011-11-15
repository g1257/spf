
/** \ingroup SPF */
/*@{*/

/*! \file Engine.h
 *
 *  SPF Engine
 *
 */
#ifndef SAVE_CONFIGS_h
#define SAVE_CONFIGS_h
#include <fstream>
#include <iostream>
#include "IoSimple.h"
#include "ProgressIndicator.h" //in PsimagLite

namespace Spf {
	
	template<typename ParametersType,typename DynVarsType>
	class SaveConfigs {
		
	public:
			
		SaveConfigs(const ParametersType& params,const DynVarsType& dynvars,size_t parallelRank)
		: params_(params),
		  dynVars_(dynvars),
		  enabled_(parallelRank==0 && params_.saveEach>0),
		  ioOut_(params_.filename+".configs",(enabled_) ? parallelRank: 1)
		{
			if (!enabled_) return;
			writeHeader();
		}
				
		~SaveConfigs()
		{
			if (!enabled_) return;
			time_t t = time(0);
			std::string s(ctime(&t));
			ioOut_<<s;
			ioOut_<<"#EOF\n";
		}
		
		void operator()(size_t iter)
		{
			if (!enabled_) return;
			if (iter>0 && iter % params_.saveEach !=0) return;
			ioOut_<<"#Iteration"<<iter<<"\n";
			ioOut_<<dynVars_;
			time_t t = time(0);
			std::string s(ctime(&t));
			ioOut_<<s;
		}
		
	private:

		
		void writeHeader()
		{
			ioOut_<<"#This is SPF v7\n";
			time_t t = time(0);
			std::string s(ctime(&t));
			ioOut_<<s;
			ioOut_<<params_;
		}

		const ParametersType& params_;
		const DynVarsType& dynVars_;
		bool enabled_;
		PsimagLite::IoSimple::Out ioOut_;
	}; // class SaveConfigs
} // namespace Spf

/*@}*/
#endif // SAVE_CONFIGS_h
