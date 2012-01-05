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
#include "TypeToString.h"

// FIXME: values_ should be a string here
namespace Spf {
	template<typename RealType,typename IoOutputType,typename ConcurrencyType>
	class  Packer {

		typedef typename ConcurrencyType::CommType CommType;

		enum {TYPE_REAL,TYPE_COMPLEX,TYPE_SIZE_T,TYPE_STRING};

	public:
		
		Packer(IoOutputType& fout,ConcurrencyType& concurrency,CommType comm)
		: fout_(fout),concurrency_(concurrency),comm_(comm),progress_("Packer",0)
		{}
		
		~Packer()
		{
			size_t nprocs = concurrency_.nprocs(comm_);
			size_t r = concurrency_.rank(comm_);
			std::vector<RealType> values(values_.size()*nprocs,0);
			for (size_t i=0;i<values_.size();i++) {
				values[r+i*nprocs] = values_[i];
			}
			values_.clear();
			
			concurrency_.reduce(values,comm_);
			
			if (!concurrency_.root(comm_)) return;
			
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

		void pack(const std::string& label,const std::string& value)
		{
			labels_.push_back(label);
			values_.push_back(atof(value.c_str()));
			types_.push_back(TYPE_STRING);
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

		IoOutputType& fout_;
		ConcurrencyType& concurrency_;
		CommType comm_;
		PsimagLite::ProgressIndicator progress_;
		std::vector<std::string> labels_;
		std::vector<RealType> values_;
		std::vector<size_t> types_;
	}; // Packer
} // namespace PsimagLite

#endif // PACKER_H

