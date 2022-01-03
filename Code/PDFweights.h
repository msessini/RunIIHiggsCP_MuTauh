/*
 * PDFweights.h
 *
 *  Created on: Sep 9, 2014
 *      Author: nehrkorn
 */

#ifndef PDFWEIGHTS_H_
#define PDFWEIGHTS_H_

#include "TString.h"
#include "LHAPDF/LHAPDF.h"

class PDFweights{
public:
	PDFweights(const char* pdfname1, const char* pdfname2){
		_pdfname1 = pdfname1;
		LHAPDF::initPDFSet(1,_pdfname1.Data());
		_pdfname2 = pdfname2; // typical PDFs: CT10nnlo, MSTW2008nlo68cl, NNPDF23_nnlo_as_0119_100
		LHAPDF::initPDFSet(2,_pdfname2.Data());
	}
	virtual ~PDFweights(){};

	int numberOfMembers(){ return LHAPDF::numberPDF(2)+1; };

	// Function to reweight event to reference set "pdfname", member "member"
	//
	// For NNPDF2x_100 one should use all 100 members and calculate the rms of the 100 values
	// of the physical observable via reweighting
	//
	// For CTEQ or MSTW, in order to assign systematics one should use the master formulae.
	// See for instance: http://cmsdoc.cern.ch/cms/PRS/gentools/www/pdfuncert/uncert.html
	double weight(int id1, int id2, double x1, double x2, double scale, unsigned int member){
		if(_pdfname2==_pdfname1) return 1.;
		LHAPDF::usePDFMember(1,0);
		double pdf1 = LHAPDF::xfx(1,x1,scale,id1)*LHAPDF::xfx(1,x2,scale,id2);
		if(pdf1<=0) return 1.;
		LHAPDF::usePDFMember(2,member);
		double pdf2 = LHAPDF::xfx(2,x1,scale,id1)*LHAPDF::xfx(2,x2,scale,id2);
		return pdf2/pdf1;
	}

private:
	TString _pdfname1;
	TString _pdfname2;
};

#endif /* PDFWEIGHTS_H_ */
