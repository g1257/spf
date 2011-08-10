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
	template<typename RealType,typename ConcurrencyType>
	class  Packer {
		enum {TYPE_REAL,TYPE_COMPLEX};
	public:
		
		Packer(size_t iter,std::ostream& fout,ConcurrencyType& concurrency)
		: iter_(iter),fout_(fout),concurrency_(concurrency),progress_("Packer",0)
		{}
		
		~Packer()
		{
			size_t nprocs = concurrency_.nprocs();
			size_t r = concurrency_.rank();
			std::vector<RealType> values(values_.size()*nprocs);
			for (size_t i=0;i<values_.size();i++) {
				values[r+i*nprocs] = values_[i];
			}
			values_.clear();
			
			concurrency_.reduce(values);
			
			if (!concurrency_.root()) return;
			
			std::string s = "iter=" + ttos(iter_);
			progress_.printline(s,fout_);
			bool flag = false;
			RealType prev = 0;
			for (size_t i=0;i<labels_.size();i++) {
				s = labels_[i] + ttos(values[i*nprocs]);
				if (types_[i] == TYPE_COMPLEX) {
					if (!flag) {
						flag = true;
						prev = values[i*nprocs];
						continue;
					} else {
						flag = false;
						std::complex<RealType> temp(values[i*nprocs],prev);
						s = labels_[i] + ttos(temp);
						prev = 0;
					}
				}
				progress_.printline(s,fout_);
			}
		}
		
		void pack(const std::string& label,const RealType& value)
		{
			labels_.push_back(label);
			values_.push_back(value);
			types_.push_back(TYPE_REAL);
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
		std::ostream& fout_;
		ConcurrencyType& concurrency_;
		PsimagLite::ProgressIndicator progress_;
		std::vector<std::string> labels_;
		std::vector<RealType> values_;
		std::vector<size_t> types_;
		
	}; // Packer
} // namespace PsimagLite

#endif // PACKER_H

