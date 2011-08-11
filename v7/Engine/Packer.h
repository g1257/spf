// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK
#ifndef PACKER_H
#define PACKER_H
#include "ProgressIndicator.h"
#include <string>
#include <iostream>
#include <fstream>

namespace Spf {
	template<typename RealType,typename IoOutputType,typename ConcurrencyType>
	class  Packer {
		enum {TYPE_REAL,TYPE_COMPLEX,TYPE_SIZE_T};
	public:
		
		Packer(IoOutputType& fout,ConcurrencyType& concurrency)
		: fout_(fout),concurrency_(concurrency),progress_("Packer",0)
		{}
		
		~Packer()
		{
			size_t nprocs = concurrency_.nprocs();
			size_t r = concurrency_.rank();
			std::vector<RealType> values(values_.size()*nprocs,0);
			for (size_t i=0;i<values_.size();i++) {
				values[r+i*nprocs] = values_[i];
			}
			values_.clear();
			
			concurrency_.reduce(values);
			
			if (!concurrency_.root()) return;
			
			bool flag = false;
			RealType prev = 0;
			for (size_t r=0;r<nprocs;r++) {
				for (size_t i=0;i<labels_.size();i++) {
					RealType val = values[r+i*nprocs];
					std::string s = labels_[i] + ttos(val);
					if (types_[i] == TYPE_COMPLEX) {
						if (!flag) {
							flag = true;
							prev = val;
							continue;
						} else {
							flag = false;
							std::complex<RealType> temp(prev,val);
							s = labels_[i] + ttos(temp);
							prev = 0;
						}
					} else if (types_[i] == TYPE_SIZE_T) {
						size_t temp = size_t(val);
						s = labels_[i] + ttos(temp);
					}
					progress_.printline(s,fout_);
				}
			}
		}
		
		void pack(const std::string& label,const RealType& value)
		{
			labels_.push_back(label);
			values_.push_back(value);
			types_.push_back(TYPE_REAL);
		}

		void pack(const std::string& label,const size_t& value)
		{
			labels_.push_back(label);
			RealType v = value;
			values_.push_back(v);
			types_.push_back(TYPE_SIZE_T);
		}

		void pack(const std::string& label,const std::complex<RealType>& value)
		{
			labels_.push_back(label);
			values_.push_back(std::real(value));
			types_.push_back(TYPE_COMPLEX);
			labels_.push_back(label);
			values_.push_back(std::imag(value));
			types_.push_back(TYPE_COMPLEX);
		}
		
	private:
		size_t iter_;
		IoOutputType& fout_;
		ConcurrencyType& concurrency_;
		PsimagLite::ProgressIndicator progress_;
		std::vector<std::string> labels_;
		std::vector<RealType> values_;
		std::vector<size_t> types_;
		
	}; // Packer
} // namespace PsimagLite

#endif // PACKER_H

