/*
 * UncertaintyValue.h
 *
 *  Created on: Jul 16, 2015
 *      Author: cherepanov
 *
 *      This class allows to store a value and an associated uncertainty
 *      in one object. Furthermore it provides operators to calculate
 *      using these objects, which propagate the errors appropriately.
 *
 *      While the value itself is templated and can be set as necessary,
 *      the uncertainty is always stored in double precision.
 *
 */

#ifndef UNCERTAINTYVALUE_H_
#define UNCERTAINTYVALUE_H_

#include <cmath>
#include <iostream>
#include <sstream>
#include <string.h>
#include "SimpleFits/FitSoftware/interface/Logger.h"

namespace UncVal{

	template <typename T>
	class UncertaintyValue {
	private:
		T val;
		double unc;

	public:
		// default constructor
		UncertaintyValue() {val = 0; unc = 0;}
		// constructor providing value and uncertainty
		UncertaintyValue(T value, double error = 0) {val = value; unc = error;}

		// obtain UncertaintyValue object with Poisson error
		static UncertaintyValue<T> Poisson(T value) {
			UncertaintyValue<T> uv(value);
			uv.poiss();
			return uv;
		}

		// destructor
		virtual ~UncertaintyValue(){
			// nothing to do
		}

		// set Poisson error
		void poiss(){
			unc = ::sqrt(val);
		}

		// getters
		T value() const {	return val;	}
		double uncer() const {	return unc;	}
		double error() const {	return uncer();}

		// arithmetic self-assignment operators
		// these should work with UncertaintyValue objects, but also with int, double, float, etc.
		template <typename TIn>
		UncertaintyValue<T>& operator+=(const TIn& right) {
			UncertaintyValue<T> r(right);
			unc = ::sqrt(unc*unc + r.error()*r.error());
			val += r.value();
			return *this;
		}
		template <typename TIn>
		UncertaintyValue<T>& operator-=(const TIn& right) {
			UncertaintyValue<T> r(right);
			unc = ::sqrt(unc*unc + r.error()*r.error());
			val -= r.value();
			return *this;
		}
		template <typename TIn>
		UncertaintyValue<T>& operator*=(const TIn& right) {
			UncertaintyValue<T> r(right);
			unc = ::sqrt(val*val*r.error()*r.error() + r.value()*r.value()*unc*unc);
			val *= r.value();
			return *this;
		}
		template <typename TIn>
		UncertaintyValue<T>& operator/=(const TIn& right) {
			UncertaintyValue<T> r(right);
			if (r.value() == 0) {
				Logger(Logger::Error) << "Division by zero." << std::endl;
				val = 0;
				unc = 0;
			}
			else {
				T lv = val;
				T rv = r.value();
				double lu = unc;
				double ru = r.uncer();
				unc = sqrt( (lu*lu)/(rv*rv) + (lv*lv*ru*ru)/(rv*rv*rv*rv));

				val /= r.value();
			}
			return *this;
		}

		// obtain numbers as string
		std::string getString() const {
			std::ostringstream stream;
			stream << "[" << value() << " +/- " << error() << "]";
			std::string s = stream.str();
			return s;
		}
	};

	// arithmetic operators
	template <typename T>
	UncertaintyValue<T> operator+(UncertaintyValue<T> left, const UncertaintyValue<T>& right) {
		left += right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator+(UncertaintyValue<T> left, const T& right) {
		left += right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator+(const T& left, UncertaintyValue<T> right) {
		right += left;
		return right;
	}
	template <typename T>
	UncertaintyValue<T> operator-(UncertaintyValue<T> left, const UncertaintyValue<T>& right) {
		left -= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator-(UncertaintyValue<T> left, const T& right) {
		left -= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator-(const T& left, const UncertaintyValue<T>& right) {
		UncertaintyValue<T> uvleft(left);
		uvleft -= right;
		return uvleft;
	}
	template <typename T>
	UncertaintyValue<T> operator*(UncertaintyValue<T> left, const UncertaintyValue<T>& right) {
		left *= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator*(UncertaintyValue<T> left, const T& right) {
		left *= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator*(const T& left, UncertaintyValue<T> right) {
		right *= left;
		return right;
	}
	template <typename T>
	UncertaintyValue<T> operator/(UncertaintyValue<T> left, const UncertaintyValue<T>& right) {
		left /= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator/(UncertaintyValue<T> left, const T& right) {
		left /= right;
		return left;
	}
	template <typename T>
	UncertaintyValue<T> operator/(const T& left, const UncertaintyValue<T>& right) {
		UncertaintyValue<T> uvleft(left);
		uvleft /= right;
		return uvleft;
	}

	// comparison operators
	template <typename T>
	bool operator==(const UncertaintyValue<T>& left, const UncertaintyValue<T>& right) {
		return left.value() == right.value() && left.uncer() == right.uncer();
	}
	template <typename T>
	bool operator!=(const UncertaintyValue<T>& left, const UncertaintyValue<T>& right) {
		return !operator==(left,right);
	}

	// output operator
	template <typename T>
	std::ostream& operator<<(std::ostream& os, const UncertaintyValue<T>& right) {
		os << right.getString();
		return os;
	}

	// math functions
	// Info: We use :: to access the global namespace
	template <typename T>
	UncertaintyValue<T> exp(const UncertaintyValue<T>& uv){
		return UncertaintyValue<T>(::exp(uv.value()), ::exp(uv.value())*uv.error() ) ;
	}
	template <typename T>
	UncertaintyValue<T> sqrt(const UncertaintyValue<T>& uv){
		return UncertaintyValue<T>(::sqrt(uv.value()), 0.5*uv.error()/::sqrt(uv.value()) ) ;
	}
	template <typename T, typename P>
	UncertaintyValue<T> pow(const UncertaintyValue<T>& uv, const P& expo){
		return UncertaintyValue<T>( ::pow(uv.value(), expo), expo*::pow(uv.value(), expo-1)*uv.error() ) ;
	}
	template <typename T>
	UncertaintyValue<T> log(const UncertaintyValue<T>& uv){
		return UncertaintyValue<T>(::log(uv.value()), uv.error()/uv.value() ) ;
	}
	template <typename T>
	UncertaintyValue<T> sin(const UncertaintyValue<T>& uv){
		return UncertaintyValue<T>(::sin(uv.value()), ::cos(uv.value())*uv.error() ) ;
	}
	template <typename T>
	UncertaintyValue<T> cos(const UncertaintyValue<T>& uv){
		return UncertaintyValue<T>(::cos(uv.value()), ::sin(uv.value())*uv.error() ) ;
	}
}

// typedefs
typedef UncVal::UncertaintyValue<double> UncDouble;
typedef UncVal::UncertaintyValue<float> UncFloat;
typedef UncVal::UncertaintyValue<int> UncInt;

#endif /* UNCERTAINTYVALUE_H_ */
