//
//  PNVmodel.cpp
//  MC2project
//
//  Created by David R Conklin on 7/18/14.
//

#include <iostream>
#include <math.h> 
#include <vector>

#include "netcdf.h"
#include "category_bgc.h" 
#include "MAPSSvegClasses.h"
#include "ScienceFcns.h"
#include "ProcessModel.h"
#include "MC2.h"
#include "MCbiogeog.h"
#include "PNVmodel.h"
#include "MAPSSbiogeographyModel.h"


PNVmodel::PNVmodel(Simulation * sP, BiogeographyInputData * inValsP)
{ // constructor   

	ProcessModelInit(sP);

	// inputs to the subregion specific routines
	prcpgrid = NC_FILL_FLOAT; // precipitation at sea level, in H2O
	foggrid = NC_FILL_FLOAT; // fog effect, each unit = 20" H2O
	cadgrid = NC_FILL_FLOAT; // cold air drainage effect, ft of elevation
	msltgrid = NC_FILL_FLOAT; // mean air temperature at sea level, deg F
	tmgrid = NC_FILL_FLOAT; // topographic moisture effect, ft of elevation
	singrid = NC_FILL_FLOAT; // aspect 
	solargrid = NC_FILL_FLOAT; // solar radiation 
	elevtmp = NC_FILL_FLOAT; // elevation, ft

	// outputs from the subregion specific routines
	vgtmp2 = NC_FILL_FLOAT; // alpine zone lower boundary, ft
	vgtmp3 = NC_FILL_FLOAT; // parkland zone lower boundary, ft
	vgtmp4 = NC_FILL_FLOAT; // mountain hemlock zone lower boundary, ft
	vgtmp5 = NC_FILL_FLOAT; // Pacific silver fir zone lower boundary, ft
	vgtmp6 = NC_FILL_FLOAT; // Sitka spruce zone upper boundary, ft
	vgtmp7 = NC_FILL_FLOAT; // subalpine fir zone lower boundary, ft
	vgtmp9 = NC_FILL_FLOAT; // Douglas-fir zone upper boundary, ft
	vgtmp10 = NC_FILL_FLOAT; // Douglas-fir zone lower boundary, ft
	vgtmp11 = NC_FILL_FLOAT; // mountain hemlock upper boundary, ft
	vgtmp12 = NC_FILL_FLOAT; // western hemlock lower boundary, ft
	vgtmp14 = NC_FILL_FLOAT; // Grand fir zone lower boundary, ft

	double elev = inValsP->elev; // elevation, m
	double pslM = 3.14066e-7; // M
	double pslN = 0.99908387; // N
	double pslP = -139015.01; // P
	double pslQ = 1922.4891; // Q
	double pslR = 561677190.; // R
	double pslS = -13825146.; // S

	foggrid = inValsP->fog; 

	tmgrid = inValsP->topomoist;

	if (!(inValsP->deltaTsl>-100.f && inValsP->deltaTsl<100.f)) 
	{
		// printf("*** PNVmodel::pnv(): deltaTsl=%f deg F, tmp_yr=%f deg C\n", inValsP->deltaTsl, inValsP->tmp_yr);
		msltgrid = NC_FILL_FLOAT;
	}
	else msltgrid = sciFn.CtoF(inValsP->tmp_yr) + inValsP->deltaTsl;

	if (!(inValsP->aspect>-100.f)) singrid = NC_FILL_FLOAT;
	else singrid = sin((inValsP->aspect - 120.f)*(PI/180.f));

	solargrid = inValsP->sw; 

	cadgrid = inValsP->cad;  
	// cad = sciFn.ft_to_m(cad_ft);  

	double tap = inValsP->ppt_yr; // total annual precipitation, mmH2O
	float psl = (float)(1./(pslM + pslN/tap + (1./(pslP + pslQ*tap))*elev + (1./(pslR + pslS*tap))*elev*elev));
	prcpgrid = sciFn.mm_to_in(psl);

	elevtmp = sciFn.m_to_ft(inValsP->elev);

}; // end of PNVmodel constructor


bool PNVmodel::runModelAndCollectData(const int year_index, const int row_index, const int col_index) 
{
	return(true);
}; // end of PNVmodel::runModelAndCollectData()


PNVcode PNVmodel::pnv(BiogeographyInputData * inValsP, ModelParamsClass * mP)
{  
	if (!(foggrid>-100.f)) return(UNKNOWNpnv); // foggrid = 0.85; // 
	if (!(tmgrid>-100.f)) return(UNKNOWNpnv); // tmgrid = 5.41; // 
	if (!(msltgrid>-100.f && msltgrid<100.f)) return(UNKNOWNpnv); // msltgrid = sciFn.CtoF(inValsP->tmp_yr + 1.f); 
	if (singrid<-1.f || singrid>1.f) return(UNKNOWNpnv); // singrid = sin((135.f - 120.f)*(PI/180.f)); // 
	if (!(solargrid>-100.f)) return(UNKNOWNpnv); // solargrid = 15048.f; // 
	if (!(cadgrid>-100.f)) return(UNKNOWNpnv); // cadgrid = 146.16; // 
	if (!(prcpgrid>=0. && prcpgrid<10000.)) return(UNKNOWNpnv);

	PNVcode pnv_code = UNKNOWNpnv;
	if (foggrid>-100.f and tmgrid>-100.f and msltgrid>-100.f && msltgrid<100.f and
			-1.f<=singrid && singrid<=1.f and solargrid>-100.f and cadgrid>-100.f and
			prcpgrid>=0.f && prcpgrid<10000.f) pnv_code = SubregionSpecificCode(inValsP);

	pS->VarDict[MATSL].save(msltgrid, YEAR_INTERVAL);
	pS->VarDict[PSL].save(prcpgrid, YEAR_INTERVAL);
	pS->VarDict[SSZ_UB].save(sciFn.ft_to_m(vgtmp6), YEAR_INTERVAL);
	pS->VarDict[PSFZ_LB].save(sciFn.ft_to_m(vgtmp5), YEAR_INTERVAL);
	pS->VarDict[SAFZ_LB].save(sciFn.ft_to_m(vgtmp7), YEAR_INTERVAL);
	pS->VarDict[PKLZ_LB].save(sciFn.ft_to_m(vgtmp3), YEAR_INTERVAL);

	return(pnv_code);
} // end of pnv()


PNVcode PNVmodel::SubregionSpecificCode(BiogeographyInputData * inValsP)
{ PNVcode pnv_code;

	switch (inValsP->subregion)
	{
		case 10010103: pnv_code = Subregion10010103(); break; // Olympic National Forest/ SW WASH
		case 10010104: pnv_code = Subregion10010104(); break; // NW Oly Peninsula
		case 10010203: pnv_code = Subregion10010203(); break; // SW WA - Willapa Hills
		case 10010204: pnv_code = Subregion10010204(); break; // Olympic NF east side   
		case 10010205: pnv_code = Subregion10010205(); break; // Mt. Hood and north Willamette NFs

		case 10020204: // 0.1% of WWC, just east of 10010206
		case 10010206: pnv_code = Subregion10010206(); break; // Gifford Pinchot N.F.

		case 10020202: // 0.2% of WNC, just east of 10010207
		case 10010207: pnv_code = Subregion10010207(); break; // Mt. Baker - Snoqualmie N.F. in WNC

		case 10010208: pnv_code = Subregion10010208(); break; // NE Oly Peninsula 
		case 10010210: pnv_code = Subregion10010210(); break; // MBS and GP National Forests, plus MRNP.  Mostly WWC
		case 10010214: pnv_code = Subregion10010214(); break; // N. Willamette Valley
		case 10010215: pnv_code = Subregion10010215(); break; // Puget lowlands

		case 10020102: // 0.1% of WNC, just east of 10020101
		case 10020110: // 0.2% of WNC, just south of 10020101 // far north Wenatchee NF, Okanogan NF
		case 10020101: pnv_code = Subregion10020101(); break; // Okanogan N.F.

		case 10020205: // A small piece in the SE corner of WWC?
		default:
			       pnv_code = UNKNOWNpnv; 
			       break;
	} // end of switch (subregion)

	return(pnv_code);
} // end of MC_BiogeographyModel::SubregionSpecificCode(PNVcode subregion)


PNVcode PNVmodel::Subregion10010104()
{ // NW corner of Olympic Peninsula
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		float _blalp1 = 100;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;
		//fog
		float _blalp12 = -250.0;
		//cad
		float _blalp11 = -0.50;

		vgtmp2 = int((((_blalp12 * foggrid + _blalp11 * cadgrid + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (0.50 * _blalp5 + _blalp6 * prcpgrid) * 
				singrid + (0.50 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		//OLY Regional
		float _blpk1 = -1000;
		float _blpk2 = 149;
		float _blpk3 = -53;
		float _blpk4 = -17;
		float _blpk5 = 250;
		float _blpk6 = -2;
		//fog
		float _blpk12 = -250.0;
		//cad
		float _blpk11 = -0.50;

		vgtmp3 = int((((_blpk12 * foggrid + _blpk11 * cadgrid + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		//OLY Regional
		float _blmh1 = -4100;
		float _blmh2 = 193;
		float _blmh3 = -181;
		float _blmh4 = -9;
		float _blmh5 = -4.0;
		float _blmh6 = 450;
		float _blmh7 = -0.8;
		float _blmh8 = -0.7;
		float _blmh9 = -0.2;
		//fog
		float _blmh12 = -300.0;
		//cad
		float _blmh11 = -0.50;

		vgtmp4 = int(((((_blmh12 * foggrid + _blmh11 * cadgrid) + 
							_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .5)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.40 * _blmh6 + _blmh7 
				  * tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid + (0.65 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		//OLY Regional
		float _blsf1 = -2900;
		float _blsf2 = 150;
		float _blsf3 = -100;
		float _blsf4 = -20;
		float _blsf5 = -0.001;
		float _blsf6 = 450;
		float _blsf7 = 7;
		float _blsf8 = 1;
		float _blsf9 = -0.35;
		//fog
		float _blsf12 = -600.0;
		//cad
		float _blsf11 = -0.55;

		vgtmp5 = int(((((_blsf12 * foggrid) + (_blsf11 * cadgrid) + _blsf1 + 
							_blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + 
				((0.40 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid + (0.65 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	{
		float _buss1 = 2200;
		float _buss2 = -115;
		float _buss3 = 40;
		float _buss4 = 250;
		float _buss5 = 1;
		float _buss6 = -350;
		float _buss7 = -20;
		float _buss8 = 500;
		float _buss9 = -500;
		// float _buss10 = -53.1; 
		//fog
		float _buss12 = 900;
		//cad
		float _buss11 = 0.70;

		vgtmp6 = int(((((_buss12 * foggrid) + (_buss11 * cadgrid)+ 
							_buss1 + _buss2 * msltgrid) + _buss3 * tmgrid) + 
					(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5)) + 
				((0.50 * _buss6 + _buss7 * tmgrid) + 
				 (_buss8 + _buss9 / tmgrid) / prcpgrid) * singrid + 
				(0.50 * _buss6 * (solargrid - 14930) / 3380.4));

		/*
		   printf("*** Subregion10010104(): foggrid=%f, cadgrid=%f, msltgrid=%f, tmgrid=%f, prcpgrid=%f, singrid=%f, solargrid=%f\n,",
		   foggrid, cadgrid, msltgrid, tmgrid, prcpgrid, singrid, solargrid);
		   printf("_buss12 * foggrid = %f\n, _buss11 * cadgrid = %f\n,_buss1 = %f\n_buss2 * msltgrid = %f\n_buss3 * tmgrid = %f\n"
		   "(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5) = %f\n((0.50 * _buss6 + _buss7 * tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) * singrid = %f\n"
		   "(0.50 * _buss6 * (solargrid - 14930) / 3380.4)) = %f\n",
		   _buss12 * foggrid, _buss11 * cadgrid, _buss1, _buss2 * msltgrid, _buss3 * tmgrid, (_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5),
		   ((0.50 * _buss6 + _buss7 * tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) * singrid, (0.50 * _buss6 * (solargrid - 14930) / 3380.4));
		   */
	}
	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 7900;
		float _blsaf2 = 100;
		float _blsaf3 = 3;
		float _blsaf4 = 9;
		float _blsaf5 = 150;
		float _blsaf6 = -30;
		float _blsaf7 = -8;
		float _blsaf8 = -0.3;
		float _blsaf9 = -75;
		//fog
		float _blsaf12 = 100.00;
		//cad
		float _blsaf11 = -0.7;

		vgtmp7 = int((( _blsaf12 * foggrid + _blsaf11 * cadgrid + 
						_blsaf1 + _blsaf9 * msltgrid + _blsaf2 * pow(tmgrid, .5)) + 
					(_blsaf3 + _blsaf4 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.50 *  _blsaf5 + _blsaf6 * tmgrid) + (_blsaf7 + _blsaf8 * 
					tmgrid) * prcpgrid) * singrid + 
				(0.50 * _blsaf5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//*******************************************************************
	{
		float _budf1 = -6900;
		float _budf2 = -10;
		float _budf3 = -.0025;
		float _budf4 = -.002;
		float _budf5 = 700;
		float _budf6 = -140;
		float _budf7 = 5000;
		float _budf8 = -194;
		float _budf9 = 150;
		//fog
		float _budf12 = -900.00;
		//cad
		float _budf11 = -0.3;

		vgtmp9 = int(_budf12 * foggrid + _budf11 * cadgrid + _budf9 * msltgrid + 
				_budf1 + _budf2 * pow(tmgrid, 1) + (_budf3 + 
					_budf4 * pow(tmgrid, 2)) * pow(prcpgrid, 3) + (0.50 * _budf5 + _budf6 * pow(tmgrid, .5) + 
					_budf7 + _budf8 * pow(tmgrid, .5) / (prcpgrid + 10)) 
				* singrid + (0.50 * _budf5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 10900;
		float _bldf2 = -222;
		float _bldf3 = 125;
		float _bldf4 = 55;
		float _bldf5 = -3;
		float _bldf6 = -900;
		float _bldf7 = 300;
		float _bldf8 = -800.0;
		float _bldf9 = -7000;
		//fog
		float _bldf12 = 1000.00;
		//cad
		float _bldf11 = 0.5;

		vgtmp10 = int((_bldf12 * foggrid + _bldf11 * cadgrid + _bldf1 + _bldf2 * msltgrid + 
					_bldf3 * pow(tmgrid, 2) + 
					(_bldf4 + _bldf5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.50 *  _bldf6 + _bldf7 * pow(tmgrid, .5)) + (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 10)) * singrid + 
				(0.50 * _bldf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 130;
		float _blgf2 = -14400;
		float _blgf3 = 406800;
		float _blgf4 = -130;
		float _blgf5 = -45;
		float _blgf6 = -9;
		float _blgf7 = 200;
		float _blgf8 = 12;
		float _blgf9 = 32;
		float _blgf10 = -7;
		//cad
		// float _blgf11 = -0.1;
		//fog
		// float _blgf12 = -100;

		vgtmp14 = int((((_blgf1 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid) + (_blgf5 + _blgf6 * tmgrid) 
					* prcpgrid) + ((0.60 * _blgf7 + _blgf8 * tmgrid) + 
						(_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
				(0.40 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary 
	//*******************************************************************
	{
		//OLY Regional
		float _blwh1 = -5600;
		float _blwh2 = 150;
		float _blwh3 = -129; 
		float _blwh4 = -40;
		float _blwh5 = -7.5;
		float _blwh6 = 750;
		float _blwh7 = -28;
		float _blwh8 = 60;
		float _blwh9 = -15;
		//CAD
		float _blwh11 = -0.8;
		//FOG
		float _blwh12 = -900.0;

		vgtmp12 = int((((_blwh12 * foggrid + _blwh11 * cadgrid + 
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * tmgrid) + (_blwh4 + 
							_blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	/*
	   printf("*** Subregion10010104(): elevtmp = %f\n", elevtmp);
	   printf("vgtmp2=%f, vgtmp3=%f, vgtmp4=%f, vgtmp5=%f, vgtmp6=%f, vgtmp7=%f, vgtmp9=%f, vgtmp10=%f, vgtmp12=%f, vgtmp14=%f\n",
	   vgtmp2, vgtmp3, vgtmp4, vgtmp5, vgtmp6, vgtmp7, vgtmp9, vgtmp10, vgtmp12, vgtmp14);
	   */

	int v14ecoregion = con ((( elevtmp < vgtmp12) and ( elevtmp > vgtmp14) and ( elevtmp < vgtmp5)), 16, 33); // GF
	int v13ecoregion = con ((( elevtmp < vgtmp4) and ( elevtmp > vgtmp12) and ( elevtmp < vgtmp5)), 19, 33); // WH
	int v11ecoregion = con ((( elevtmp < vgtmp12 and ( elevtmp < vgtmp5)and ( elevtmp < vgtmp10)  )), 4, 33); // GL
	int v10ecoregion = con (((elevtmp < vgtmp6)), 9, 33); // SS
	// int v9ecoregion = con (elevtmp < vgtmp5, 19, 33); // WH
	int v8ecoregion = con (((elevtmp < vgtmp9) and (elevtmp > vgtmp10) and (elevtmp < vgtmp3)), 14, 33); // DF
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp3)), 22,  33); // PSF
	int v6ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp2)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp7) and (elevtmp < vgtmp2)), 1, 33); // SAF
	int v4ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

	int vz_ecoregion2 = v4ecoregion;
	if (v5ecoregion<vz_ecoregion2) vz_ecoregion2 = v5ecoregion;
	if (v6ecoregion<vz_ecoregion2) vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v8ecoregion<vz_ecoregion2) vz_ecoregion2 = v8ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;
	if (v13ecoregion<vz_ecoregion2) vz_ecoregion2 = v13ecoregion;
	if (v14ecoregion<vz_ecoregion2) vz_ecoregion2 = v14ecoregion;


	//---Renumber saf to match numbers in table and shadesets----------
	int vz_ecoregion = con(vz_ecoregion2 == 1, 25, vz_ecoregion2);

	return((PNVcode)vz_ecoregion);

} // end of Subregion10010104()




PNVcode PNVmodel::Subregion10010203()
{ // SW Washington-Willapa Hills

	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		//-----Alpine----------
		int _blalp1 = 200;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;
		vgtmp2 = int((((_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (_blalp5 + _blalp6 * prcpgrid) * singrid);
	}
	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		//-----Parkland-------
		float _blpk1 = -1400;
		float _blpk2 = 150; 
		float _blpk3 = -50;
		float _blpk4 = -10;
		float _blpk5 = 180;
		float _blpk6 = 2;

		vgtmp3 = int((((_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
	}
	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	//----Lower MH-----
	{
		float _blmh1 = -3720;
		float _blmh2 = 186;
		float _blmh3 = -234;
		float _blmh4 = -10;
		float _blmh5 = -3.1;
		float _blmh6 = 450;
		float _blmh7 = 2;
		float _blmh8 = -.5;
		float _blmh9 = -.2;

		vgtmp4 = int((((_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .5)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((0.50 * _blmh6 + _blmh7 
							* tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid + (0.50 * _blmh6 * (solargrid - 14930) / 3380.4));
	}
	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	//----Lower PSF-----
	{
		float _blsf1 = -5900;
		float _blsf2 = 245;
		float _blsf3 = -420;
		float _blsf4 = -26;
		float _blsf5 = .22;
		float _blsf6 = 400;
		float _blsf7 = 6;
		float _blsf8 = .4;
		float _blsf9 = -.35;

		vgtmp5 = int((((_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + ((0.50 * _blsf6 + 
							_blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid + (0.50 * _blsf6 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary  from 10103  02-18-05
	//*******************************************************************
	//----Upper SS--from 10103--021805---
	//float _buss1 10200
	{
		float _buss1 = 10000;
		float _buss2 = -303;
		float _buss3 = 820;
		float _buss4 = 210;
		float _buss5 = -0.75;
		float _buss6 = -10;
		float _buss7 = -25;
		float _buss8 = 500;
		float _buss9 = -500;
		float _buss10 = -52.9;
		float _buss12 = 100;
		float _buss11 = 1.0;

		vgtmp6 = int((((_buss12 * foggrid + _buss11 * cadgrid + _buss1 + _buss2 * msltgrid) + _buss3 * tmgrid + 
						_buss10 * pow(tmgrid, 2)) + (_buss4 + _buss5 * pow(tmgrid, 1.8)) * pow(prcpgrid, .5)) + 
				((0.50 * _buss6 + _buss7 * 
				  tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) 
				* singrid + (0.50 * _buss6 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	{
		float _blwh1 = -8500;
		float _blwh2 = 290;
		float _blwh3 = -220; 
		float _blwh4 = -80;
		float _blwh5 = -8;
		float _blwh6 = 550;
		float _blwh7 = 30;
		float _blwh8 = 70;
		float _blwh9 = -14;

		//cad
		float _blwh11 = -0.6;
		//fog
		float _blwh12 = -1200;


		vgtmp7 = int((((_blwh12 * foggrid + _blwh11 * cadgrid +
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * 
						tmgrid) + (_blwh4 + _blwh5 * tmgrid) * prcpgrid) + 
				((0.50 * _blwh6 + _blwh7 * tmgrid) + 
				 (_blwh8 + _blwh9 * tmgrid) * prcpgrid) 
				* singrid + (0.50 * _blwh6 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//*******************************************************************
	////----Upper DF---added 2-18-05--
	////float _budf1 255600
	//float _budf1 255800
	{
		float _budf1 = 254300;
		float _budf2 = -9000;
		float _budf3 = 81;
		float _budf4 = -120;
		float _budf5 = -90;
		float _budf6 = -9;
		float _budf7 = 1000;
		float _budf8 = -80;
		float _budf9 = 20000;
		float _budf10 = -200;
		//cad
		float _budf11 = -0.2;
		//fog
		float _budf12 = -1000;

		vgtmp9 = int((((_budf11 * cadgrid + _budf12 * foggrid + 
							_budf1 + _budf2 * msltgrid + _budf3 * pow(msltgrid, 2)) 
						+ _budf4 * tmgrid) + 
					(_budf5 + _budf6 * tmgrid) * prcpgrid) + ((_budf7 + _budf8 * 
							pow(tmgrid, .5)) + (_budf9 + _budf10 * tmgrid) / (prcpgrid + 10)) 
				* singrid);
	}
	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	//----Lower DF-----
	{
		float _bldf1 = 700000;
		float _bldf2 = -33900;
		float _bldf3 = 539.1;
		float _bldf4 = -11;
		float _bldf5 = -9;
		float _bldf6 = -1300;
		float _bldf7 = 1100;
		float _bldf8 = -10000;
		float _bldf9 = 100;
		//CAD
		float _bldf11 = -0.5;
		//FOG
		float _bldf12 = 300.0;
		float _bldf14 = -2.8;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf6 + _bldf7 * pow(tmgrid, .5)) + 
				 (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 10)) 
				* singrid + (0.60 * _bldf6 * (solargrid - 14930) / 3380.4));

	}
	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	//----Lower GF------added 2-18-05
	{
		float _blgf1 = 1168600;
		float _blgf2 = -43160;
		float _blgf3 = 395.5;
		float _blgf4 = -400;
		float _blgf5 = 23000;
		float _blgf6 = -500;
		float _blgf7 = 400;
		float _blgf8 = 5;
		float _blgf9 = 25;
		float _blgf10 = -5;
		float _blgf11 = -0.5;
		// float _blgf12 = -200;
		float _blgf12 = 200.0;
		float _blgf15 = 5;

		vgtmp12 = int(((_blgf11 * cadgrid + _blgf12 * foggrid + 
						_blgf3 * pow(msltgrid, 2) + _blgf2 * msltgrid +_blgf1) +  
					+ _blgf4 * tmgrid + _blgf5 + _blgf6 * prcpgrid) + _blgf15 * pow(prcpgrid, 2)) +
			((_blgf7 + _blgf8 * tmgrid) + 
			 (_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid;

	}
	//*******************************************************************


	//=====================================================================
	//vgtmp2=lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=ubss; vgtmp7=lbwh;
	//[vgtmp8=ubsaf]; vgtmp9=ubdf; vgtmp10=lbdf; [vgtmp11=ubmh]; vgtmp12=lbgf
	//=====================================================================

	int v11ecoregion = con (((elevtmp < vgtmp7) and (elevtmp < vgtmp10) and (elevtmp < vgtmp12)), 4, 22); // GL lbwh, lbgf, lbdf
	int v10ecoregion = con (elevtmp < vgtmp6, 9, 23); // SS
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 23); // WH
	int v8ecoregion = con (((elevtmp < vgtmp7) and (elevtmp < vgtmp12) /* and (elevtmp > vgtmp10) */), 14, 23); // DF lbwh, lbfg
	int v7ecoregion = con (elevtmp < vgtmp4, 22, 23); // SF
	int v6ecoregion = con (elevtmp < vgtmp3, 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp7) and (elevtmp > vgtmp12)), 16, 23); // GF lbwh, lbgf

	int vz_ecoregion2 = v5ecoregion;
	if (v6ecoregion<vz_ecoregion2) vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v8ecoregion<vz_ecoregion2) vz_ecoregion2 = v8ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;

	return((PNVcode)vz_ecoregion2);

} // end of int PNVmodel::Subregion10010203()


PNVcode PNVmodel::Subregion10010204()
{ // Olympic N.F. east side
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	//OLY Regional
	//-----Alpine----------
	{
		int _blalp1 = 100;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;
		//fog
		float _blalp12 = -250.0;
		//cad
		float _blalp11 = -0.50;

		vgtmp2 = int((((_blalp12 * foggrid + _blalp11 * cadgrid + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (0.50 * _blalp5 + _blalp6 * prcpgrid) * 
				singrid + (0.50 * _blalp5 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	//OLY Regional
	//-----Parkland-------
	{
		float _blpk1 = -1100;
		float _blpk2 = 149; 
		float _blpk3 = -53;
		float _blpk4 = -17;
		float _blpk5 = 250;
		float _blpk6 = -2;
		//fog
		float _blpk12 = -250.0;
		//cad
		float _blpk11 = -0.50;

		vgtmp3 = int((((_blpk12 * foggrid + _blpk11 * cadgrid + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
	}
	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	//OLY Regional
	//----Lower MH-----
	{
		float _blmh1 = -3400;
		float _blmh2 = 218;
		float _blmh3 = -395;
		float _blmh4 = -32;
		float _blmh5 = 0.15;
		float _blmh6 = 400;
		float _blmh7 = 2.0;
		float _blmh8 = -0.5;
		float _blmh9 = -0.2;
		//fog
		float _blmh12 = -300.0;
		//cad
		float _blmh11 = -0.55;

		vgtmp4 = int(((( (_blmh12 * foggrid + _blmh11 * cadgrid) + _blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, 1)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, 2)) * prcpgrid) + ((0.40 * _blmh6 + _blmh7 * tmgrid) + 
						(_blmh8 + _blmh9 * tmgrid) * prcpgrid) * singrid + (0.65 * _blmh6 * (solargrid - 14930) / 3380.4));
	}
	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -3900;
		float _blsf2 = 212;
		float _blsf3 = -330;
		float _blsf4 = -41;
		float _blsf5 = 0.15;
		float _blsf6 = 400;
		float _blsf7 = 6;
		float _blsf8 = 0.4;
		float _blsf9 = -.35;
		//fog
		float _blsf12 = -300.0;
		//cad
		float _blsf11 = -0.55;

		vgtmp5 = int(((((_blsf12 * foggrid) + (_blsf11 * cadgrid) + _blsf1 + 
							_blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + 
				((0.40 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid + (0.65 * _blsf6 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	{
		float _buss1 = 10300;
		float _buss2 = -278;
		float _buss3 = 66;
		float _buss4 = -20;
		float _buss5 = 40;
		float _buss6 = -70;
		float _buss7 = -20;
		float _buss8 = 500;
		float _buss9 = -500;
		//fog
		float _buss12 = 100.00;
		//cad
		float _buss11 = 0.7;

		vgtmp6 = int(((((_buss12 * foggrid) + (_buss11 * cadgrid)+ _buss1 + _buss2 * 
							msltgrid) + _buss3 * tmgrid) + 
					(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5)) + ((0.50 * _buss6 + _buss7 * 
						tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) 
				* singrid + (0.50 * _buss6 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 8100;
		float _blsaf2 = 110;
		float _blsaf3 = 3;
		float _blsaf4 = 9;
		float _blsaf5 = 150;
		float _blsaf6 = -30;
		float _blsaf7 = -8;
		float _blsaf8 = -0.3;
		float _blsaf9 = -75;
		//fog
		float _blsaf12 = 350.00;
		//cad
		float _blsaf11 = -0.7;

		vgtmp7 = int((( _blsaf12 * foggrid + _blsaf11 * cadgrid + 
						_blsaf1 + _blsaf9 * msltgrid + _blsaf2 * pow(tmgrid, .5)) + 
					(_blsaf3 + _blsaf4 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.50 *  _blsaf5 + _blsaf6 * tmgrid) + (_blsaf7 + _blsaf8 * 
					pow(tmgrid, .5)) * prcpgrid) * singrid + (0.50 * _blsaf5 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//*******************************************************************
	//OLY Regional
	//----Upper DF/LBWHZ
	{
		float _budf1 = -7000;
		float _budf2 = -100;
		float _budf3 = 0.020;
		float _budf4 = -.001;
		float _budf5 = 850;
		float _budf6 = -100;
		float _budf7 = 10000;
		float _budf8 = -140;
		float _budf9 = 200;
		//CAD
		float _budf11 = 0.9;
		//FOG
		float _budf12 = -100.0;

		vgtmp9 = int(( _budf12 * foggrid + _budf11 * cadgrid + _budf9 * msltgrid + 
					_budf1 + _budf2 * pow(tmgrid, 2) + 
					(_budf3 + _budf4 * pow(tmgrid, 2)) * pow(prcpgrid, 3)) + 
				((0.50 * _budf5 + _budf6 * 
				  pow(tmgrid, .5)) + (_budf7 + _budf8 * tmgrid) / (prcpgrid + 10))  
				* singrid + (0.50 * _budf5 * (solargrid - 14930) / 3380.4));
	}
	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		//----Lower DF--new as of 071414---from 10204-072714- and then from 10206, rev 073014
		float _bldf1 = 696500;
		float _bldf2 = -33920;
		float _bldf3 = 539.1;
		float _bldf4 = -10;
		float _bldf5 = -7;
		float _bldf6 = -1300;
		float _bldf7 = 1100;
		float _bldf8 = -10000;
		float _bldf9 = 100;
		//CAD
		float _bldf11 = -0.5;
		//FOG
		float _bldf12 = 300.0;
		float _bldf14 = -2.8;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf6 + _bldf7 * pow(tmgrid, .5)) + 
				 (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 10)) 
				* singrid + (0.60 * _bldf6 * (solargrid - 14930) / 3380.4));

	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	//OLY Regional

	//----Upper MH---
	{float _bumh1 = 3000;
		float _bumh2 = -25;
		float _bumh3 = -80000;
		float _bumh4 = 80;
		float _bumh5 = -500;
		float _bumh6 = 7;
		float _bumh8 = 400;
		float _bumh9 = -1.2;
		//CAD
		float _bumh11 = -0.15;
		//FOG
		float _bumh12 = -100.0;

		vgtmp11 = int(_bumh12 * foggrid + _bumh11 * cadgrid + _bumh8 * msltgrid + 
				_bumh9 * pow(msltgrid, 2.35) + _bumh1 + _bumh2 * tmgrid + (_bumh3 + _bumh4 * 
					pow(tmgrid, 3)) / prcpgrid + (0.50 * _bumh5 + _bumh6 * pow(tmgrid, 2)) * singrid + 
				0.50 * _bumh5 * (solargrid - 14930) / 3380.4);
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	//----Lower GF------added 7-9-14----
	{
		float _blgf1 = 364.1;
		float _blgf2 = -40175.5;
		float _blgf3 = 1123300;
		float _blgf4 = -4150;
		float _blgf5 = 1400;
		float _blgf6 = -132;
		float _blgf7 = 900;
		float _blgf8 = -60;
		float _blgf9 = -60;
		float _blgf10 = 5;
		//cad
		float _blgf11 = -0.2;
		//fog
		float _blgf12 = -800;
		float _blgf13 = 300;
		float _blgf14 = 1.6;

		vgtmp14 = int(((((_blgf12 * foggrid) + (_blgf11 * cadgrid) + 
							_blgf1 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid + _blgf13 * pow(tmgrid, 2) ) + 
					(_blgf5 + _blgf6 * prcpgrid + _blgf14 * pow(prcpgrid, 2)) 
					+ ((0.60 * _blgf7 + _blgf8 * tmgrid) +
						(_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
					0.40 * _blgf7 * (solargrid - 14930) / 3380.4));

	}

	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	//OLY Regional
	{
		float _blwh1 = -5700;
		float _blwh2 = 196;
		float _blwh3 = -180; 
		float _blwh4 = -5;
		float _blwh5 = -15;
		float _blwh6 = 600;
		float _blwh7 = -8;
		float _blwh8 = 30;
		float _blwh9 = -5;
		//CAD
		float _blwh11 = -2.5;
		//FOG
		float _blwh12 = -900.0;

		vgtmp12 = int(_blwh12 * foggrid + _blwh11 * cadgrid + 
				_blwh1 + _blwh2 * msltgrid + _blwh3 * tmgrid + (_blwh4 + 
					_blwh5 * tmgrid) * prcpgrid + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));

	}
	//=====================================================================
	//vgtmp2 = lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=ubss; vgtmp7=lbsafz; 
	//vgtmp8=ubsafz; vgtmp9=ubdfz; vgtmp10=lbdfz; vgtmp11=ubmhz; vgtmp12=lbwhz; vgtmp14=lbgfz; 

	//elev and psl version 070914

	int v14ecoregion = con ((( elevtmp < vgtmp12) and ( elevtmp < vgtmp7) and ( elevtmp < vgtmp5)), 16, 33); // GF
	int v13ecoregion = con ((( elevtmp < vgtmp12) and ( elevtmp < vgtmp14)), 14, 33); // DF
	// int v11ecoregion = con ((( elevtmp < vgtmp10 and ( elevtmp < vgtmp14) and ( elevtmp < vgtmp12)  )), 4, 33); // GL
	int v10ecoregion = con (((elevtmp < vgtmp6)), 9, 33); // SS
	int v9ecoregion = con ((( elevtmp < vgtmp4) and ( elevtmp < vgtmp5) and ( elevtmp < vgtmp7)), 19, 33); // WH
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp3)), 22,  33); // PSF
	int v6ecoregion = con (((elevtmp < vgtmp3)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp7)), 1, 33); // SAF
	int v4ecoregion = con ((elevtmp < vgtmp2) and (elevtmp > vgtmp3), 32, 33); // PK

	int vz_ecoregion2 = v4ecoregion;
	if (v5ecoregion<vz_ecoregion2) vz_ecoregion2 = v5ecoregion;
	if (v6ecoregion<vz_ecoregion2) vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	//if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;
	if (v13ecoregion<vz_ecoregion2) vz_ecoregion2 = v13ecoregion;
	if (v14ecoregion<vz_ecoregion2) vz_ecoregion2 = v14ecoregion;

	//---Renumber saf to match numbers in table and shadesets----------
	int vz_ecoregion = con(vz_ecoregion2 == 1, 25, vz_ecoregion2);

	return((PNVcode)vz_ecoregion);

} // end of int PNVmodel::Subregion10010204()


PNVcode PNVmodel::Subregion10010207()
{ // Mt. Baker - Snoqualmie N.F. in WNC
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	// create the elevation boundary grid
	{
		//-----Alpine----------
		float _blalp1 = 1300;
		float _blalp2 = 120;
		float _blalp3 = -40;
		float _blalp4 = -11;
		float _blalp5 = 200;
		float _blalp6 = 0.5;
		//cad
		float _blalp11 = -.11;
		//fog 
		float _blalp12 = -300;

		vgtmp2 = int((((_blalp12 * foggrid + _blalp11 * cadgrid + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (0.40 * _blalp5 + _blalp6 * prcpgrid) * singrid + 
				(0.60 * _blalp5 * (solargrid - 14930) / 3380.4));

	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		//-----Parkland-------
		float _blpk1 = -1100;
		float _blpk2 = 157;
		float _blpk3 = -75;
		float _blpk4 = -22;
		float _blpk5 = 450;
		float _blpk6 = -2.0;
		//cad
		float _blpk11 = -1.0;
		//fog with PSL not APSL
		float _blpk12 = -400;

		vgtmp3 = int((((_blpk12 * foggrid + _blpk11 * cadgrid + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + 
				(0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));

	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		float _blmh1 = 300;
		float _blmh2 = 128;
		float _blmh3 = -210;
		float _blmh4 = -13;
		float _blmh5 = -2.4;
		float _blmh6 = 350;
		float _blmh7 = 15;
		float _blmh8 = -0.55;
		float _blmh9 = -0.3;
		//cad
		float _blmh11 = -0.45;
		//fog
		float _blmh12 = -400;

		vgtmp4 = int(((( _blmh12 * foggrid + _blmh11 * cadgrid + 
							_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, 1)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.40 * _blmh6 + _blmh7 * tmgrid) + 
				 (_blmh8 + _blmh9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -1000;
		float _blsf2 = 133;
		float _blsf3 = -250;
		float _blsf4 = -25;
		float _blsf5 = 0.5;
		float _blsf6 = 500;
		float _blsf7 = 10;
		float _blsf8 = -0.5;
		float _blsf9 = -.25;
		//fog 
		float _blsf12 = -500;
		//cad
		float _blsf11 = -1.2;

		vgtmp5 = int(((((_blsf12 * foggrid) + (_blsf11 * cadgrid) + 
							_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.40 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		//----Lower SAF-------
		//after graph anal of 071811, rev and imported from 20202-080814
		float _blsaf1 = 7400;
		float _blsaf2 = -49;
		float _blsaf3 = -30;
		float _blsaf4 = 8;
		float _blsaf5 = 5;
		float _blsaf6 = -20;
		float _blsaf7 = 60;
		float _blsaf8 = -10;
		float _blsaf9 = -0.2;
		//cad
		float _blsaf11 = -0.4;
		//fog
		float _blsaf12 = 60;

		vgtmp7 = int((((_blsaf11 * cadgrid + _blsaf12 * foggrid + 
							_blsaf1 + _blsaf2 * msltgrid) + (_blsaf3 * pow(tmgrid, .5))) + 
					(_blsaf4 + _blsaf5 * pow(tmgrid, .5)) * prcpgrid) + 
				(0.40 * _blsaf6 + _blsaf7 * tmgrid + (_blsaf8 + _blsaf9 * 
								      tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blsaf6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//********************************************************************
	{
		float _blwh1 = -9000;
		float _blwh2 = 270;
		float _blwh3 = -150;
		float _blwh4 = -50;
		float _blwh5 = -8;
		float _blwh6 = 500;
		float _blwh7 = 25;
		float _blwh8 = 60;
		float _blwh9 = -20;
		//CAD
		float _blwh11 = -0.6;
		//FOG
		float _blwh12 = -1200.0;

		vgtmp9 = int( _blwh12 * foggrid + _blwh11 * cadgrid + 
				_blwh1 + _blwh2 * msltgrid + _blwh3 * tmgrid + 
				(_blwh4 + _blwh5 * tmgrid) * prcpgrid + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*********************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 699000;
		float _bldf2 = -33900;
		float _bldf3 = 539.1;
		float _bldf4 = -12;
		float _bldf5 = -8;
		float _bldf7 = -1300;
		float _bldf8 = 700;
		float _bldf9 = -10000;
		float _bldf10 = 200;
		//CAD
		float _bldf11 = -0.5;
		//FOG
		float _bldf12 = 300.0;
		float _bldf14 = -2.8;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf7 + _bldf8 * pow(tmgrid, .5)) + 
				 (_bldf9 + _bldf10 * tmgrid) / (prcpgrid + 10)) 
				* singrid + (0.60 * _bldf7 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		float _bumh1 = 6500;
		float _bumh3 = -75;
		float _bumh4 = -65000;
		float _bumh5 = 70;
		float _bumh6 = -800;
		float _bumh7 = 12.0;
		float _bumh2 = 10;
		//cad
		float _bumh11 = 1.3;
		//fog
		float _bumh12 = 600;

		vgtmp11 = int(((_bumh11 * cadgrid + _bumh12 * foggrid + 
						_bumh1 + _bumh2 * msltgrid + _bumh3 * tmgrid) + ( _bumh4 + 
							_bumh5 * pow(tmgrid, 3)) / prcpgrid) + ((_bumh6 + _bumh7 
							* pow(tmgrid, 2))) * singrid);
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 364.12;
		float _blgf2 = -40175.5;
		float _blgf3 = 1114500;
		float _blgf4 = -700;
		float _blgf5 = -125;
		float _blgf6 = 2;
		float _blgf7 = 450;
		float _blgf8 = 10;
		float _blgf9 = 25;
		float _blgf10 = -6;
		//cad
		float _blgf11 = -0.4;
		//fog
		float _blgf12 = -500;
		float _blgf14 = -100;
		float _blgf15 = 3;

		vgtmp12 = int(((((_blgf12 * foggrid) + (_blgf11 * cadgrid) + 
							_blgf1 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid) + (_blgf5 + _blgf6 * pow(tmgrid, 2)) 
					* prcpgrid) + _blgf14 + _blgf15 * pow(msltgrid, 2) + 
				((0.5 * _blgf7 + _blgf8 * tmgrid) + 
				 (_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
				(0.50 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//===================================ecoregion 10010207==================================

	//vgtmp2=lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=ubss; vgtmp7=lbsaf;
	//vgtmp9=ubdf/lbwh; vgtmp10=lbdf; vgtmp11=ubmh; 

	// v10%ecoregion% = con (((elevtmp < vgtmp5) and (elevtmp < vgtmp10) and (elevtmp < vgtmp9)), 4, 33)  /* GL
	int v9ecoregion = con (((elevtmp < vgtmp5) and (elevtmp < vgtmp7)), 19, 33);  // WH
	int v8ecoregion = con (((elevtmp < vgtmp9) and (elevtmp > vgtmp10)), 14, 33); // DF
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp11) and (elevtmp < vgtmp3)), 22, 33); // SF
	int v6ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp > vgtmp11) and (elevtmp < vgtmp3) and (elevtmp > vgtmp7)), 25, 33); // SAF
	int v4ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

	int vz_ecoregion = v4ecoregion;
	if (v5ecoregion<vz_ecoregion) vz_ecoregion = v5ecoregion; 
	if (v6ecoregion<vz_ecoregion) vz_ecoregion = v6ecoregion; 
	if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion; 
	if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion; 
	if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion; 

	return((PNVcode)vz_ecoregion);

} // end of Subregion10010207()


PNVcode PNVmodel::Subregion10010208()
{ // NE corner of Olympic Peninsula
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		int _blalp1 = 100;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 220;
		float _blalp6 = -2.1;
		//CAD
		float _blalp11 = -0.25;
		//FOG
		float _blalp12 = -300.0;

		vgtmp2 = int(((( (   _blalp12 * foggrid) + (_blalp11 * cadgrid) + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (0.50 * _blalp5 + _blalp6 * prcpgrid) * 
				singrid + (0.50 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -1200;
		float _blpk2 = 148 ;
		float _blpk3 = -53;
		float _blpk4 = -17;
		float _blpk5 = 150;
		float _blpk6 = -2;
		//CAD
		float _blpk11 = -0.15;
		//FOG
		float _blpk12 = -300.0;

		vgtmp3 = int(((((   _blpk12 * foggrid) + (_blpk11 * cadgrid) + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		float _blmh1 = -5450;
		float _blmh2 = 214;
		float _blmh3 = -233;
		float _blmh4 = -8;
		float _blmh5 = -3.0;
		float _blmh6 = 350;
		float _blmh7 = 2;
		float _blmh8 = -.5;
		float _blmh9 = -.2;
		//CAD
		float _blmh11 = -0.45;
		//FOG
		float _blmh12 = -700.0;

		vgtmp4 = int(((((   _blmh12 * foggrid) + (_blmh11 * cadgrid) + 
							_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .5)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((0.40 * _blmh6 + _blmh7 
							* tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid + (0.65 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -6150;
		float _blsf2 = 245;
		float _blsf3 = -435;
		float _blsf4 = -26;
		float _blsf5 = .23;
		float _blsf6 = 400;
		float _blsf7 = 6;
		float _blsf8 = .4;
		float _blsf9 = -.35;
		//FOG
		float _blsf12 = -800;
		////CAD
		float _blsf11 = -0.6;

		vgtmp5 = int(((( (_blsf12 * foggrid) + (_blsf11 * cadgrid) + _blsf1 + 
							_blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + 
				((0.40 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid + (0.65 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	{
		float _buss1 = 10300;
		float _buss2 = -278;
		float _buss3 = 66;
		float _buss4 = -20;
		float _buss5 = 40;
		float _buss6 = -70;
		float _buss7 = -20;
		float _buss8 = 500;
		float _buss9 = -500;
		//FOG
		float _buss12 = 300.0;
		//CAD
		float _buss11 = 0.7;

		vgtmp6 = int(((((_buss12 * foggrid) + (_buss11 * cadgrid) + 
							_buss1 + _buss2 * 
							msltgrid) + _buss3 * tmgrid) + 
					(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5)) + ((0.50 * 
						_buss6 + _buss7 * 
						tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) 
				* singrid + (0.50 * _buss6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 7500;
		float _blsaf2 = 100;
		float _blsaf3 = 3;
		float _blsaf4 = 12;
		float _blsaf5 = 150;
		float _blsaf6 = -30;
		float _blsaf7 = -8;
		float _blsaf8 = -0.3;
		float _blsaf9 = -75;
		//CAD
		float _blsaf11 = -0.45;
		//FOG
		float _blsaf12 = 200.0;

		vgtmp7 = int((((_blsaf12 * foggrid) + (_blsaf11 * cadgrid) + 
						_blsaf1 + _blsaf9 * msltgrid + _blsaf2 * pow(tmgrid, .5)) + 
					(_blsaf3 + _blsaf4 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.50 *  _blsaf5 + _blsaf6 * tmgrid) + (_blsaf7 + _blsaf8 * 
					pow(tmgrid, .5)) * prcpgrid) * singrid + (0.50 * _blsaf5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary //LBWHZ
	//*******************************************************************
	{
		float _budf1 = -5700;
		float _budf2 = -100;
		float _budf3 = -0.020;
		float _budf4 = -.0010;
		float _budf5 = 950;
		float _budf6 = -140;
		float _budf7 = 10000;
		float _budf8 = -140;
		float _budf9 = 200;
		//CAD
		float _budf11 = 0.9;
		//FOG
		float _budf12 = -100.0;

		vgtmp9 = int((_budf12 * foggrid + _budf11 * cadgrid + _budf9 * msltgrid + 
					_budf1 + _budf2 * pow(tmgrid, 2) + (_budf3 + 
						_budf4 * pow(tmgrid, 2)) * pow(prcpgrid, 3) + 
					((0.50 * _budf5 + _budf6 * pow(tmgrid, .5)) + 
					 (_budf7 + _budf8 * tmgrid) / (prcpgrid + 10)) * singrid + 
					(0.50 * _budf5 * (solargrid - 14930) / 3380.4)));
	}

	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 696500;
		float _bldf2 = -33920;
		float _bldf3 = 539.1;
		float _bldf4 = -10;
		float _bldf5 = -7;
		float _bldf7 = -1300;
		float _bldf8 = 1100;
		float _bldf9 = -10000;
		float _bldf10 = 100;
		//CAD
		float _bldf11 = -0.5;
		//FOG
		float _bldf12 = -400.0;
		float _bldf14 = -2.778;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf7 + _bldf8 * pow(tmgrid, .5)) + 
				 (_bldf9 + _bldf10 * tmgrid) / (prcpgrid + 10))
				* singrid + (0.60 * _bldf7 * (solargrid - 14930) / 3380.4));

	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		float _bumh1 = 3000;
		float _bumh2 = -25;
		float _bumh3 = -80000;
		float _bumh4 = 80;
		float _bumh5 = -500;
		float _bumh6 = 7;
		float _bumh8 = 400;
		float _bumh9 = -1.2;
		//CAD
		float _bumh11 = -0.15;
		//FOG
		float _bumh12 = -100.0;

		vgtmp11 = int(((( _bumh12 * foggrid) + (_bumh11 * cadgrid) + _bumh8 * msltgrid 
						+ _bumh9 * pow(msltgrid, 2.35) + _bumh1 + _bumh2 * tmgrid) + ( _bumh3 + 
						_bumh4 * pow(tmgrid, 3)) / prcpgrid) + ((0.50 * _bumh5 + _bumh6 
						* pow(tmgrid, 2))) * singrid + (0.50 * _bumh5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir      10208
	//*******************************************************************
	{
		float _blwh1 = -10000;
		float _blwh2 = 290;
		float _blwh3 = -180; 
		float _blwh4 = -55;
		float _blwh5 = -8;
		float _blwh6 = 500;
		float _blwh7 = 30;
		float _blwh8 = 20;
		float _blwh9 = -5;
		//CAD
		float _blwh11 = -0.6;
		//FOG
		float _blwh12 = -1500.0;

		vgtmp12 = int((((_blwh12 * foggrid + _blwh11 * cadgrid + 
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * tmgrid) + (_blwh4 + 
							_blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 364.10;
		float _blgf2 = -40175.5;
		float _blgf3 = 1123900;
		float _blgf4 = -4150;
		float _blgf5 = 1400;
		float _blgf6 = -132;
		float _blgf7 = 900;
		float _blgf8 = -60;
		float _blgf9 = -60;
		float _blgf10 = -5;
		//cad
		float _blgf11 = -0.2;
		//fog
		float _blgf12 = -800;
		float _blgf13 = 300;
		float _blgf14 = 1.6;

		vgtmp14 = int(((((_blgf12 * foggrid) + (_blgf11 * cadgrid) + 
							_blgf1 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid + _blgf13 * pow(tmgrid, 2) ) + 
					(_blgf5 + _blgf6 * prcpgrid + _blgf14 * pow(prcpgrid, 2)) 
					+ ((0.60 * _blgf7 + _blgf8 * tmgrid) + 
						(_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
					0.40 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//=====================================================================
	//vgtmp2 = lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=ubss; vgtmp7=lbsafz; 
	//vgtmp8=ubsafz;vgtmp9=ubdfz; vgtmp10=lbdfz; vgtmp11=ubmhz; vgtmp12=lbwhz; vgtmp13=lb2dfz; vgtmp14=lbgfz

	int v14ecoregion = con ((( elevtmp < vgtmp12) and ( elevtmp < vgtmp7) and ( elevtmp < vgtmp5)), 16, 33); // GF
	// int v11ecoregion = con ((( elevtmp < vgtmp10 and ( elevtmp < vgtmp14) and ( elevtmp < vgtmp12) )), 4, 33); // GL
	int v10ecoregion = con ((( elevtmp < vgtmp6)), 9, 33); // SS
	int v9ecoregion = con ((( elevtmp < vgtmp4) and ( elevtmp < vgtmp5) and ( elevtmp < vgtmp7)), 19, 33); // WH
	int v13ecoregion = con ((( elevtmp < vgtmp12) and ( elevtmp < vgtmp14)), 14, 33); // DF
	int v7ecoregion = con ((( elevtmp < vgtmp4) and ( elevtmp < vgtmp3)), 22,  33); // PSF
	int v6ecoregion = con ((( elevtmp < vgtmp3)), 23, 33); // MH
	int v5ecoregion = con (( ( elevtmp < vgtmp3) and ( elevtmp > vgtmp7)), 1, 33); // SAF
	int v4ecoregion = con (( elevtmp < vgtmp2 and ( elevtmp > vgtmp3)), 32, 33); // PK

	int vz_ecoregion2 = v4ecoregion;
	if (v5ecoregion<vz_ecoregion2) vz_ecoregion2 = v5ecoregion;
	if (v6ecoregion<vz_ecoregion2) vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	//if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;
	if (v13ecoregion<vz_ecoregion2) vz_ecoregion2 = v13ecoregion;
	if (v14ecoregion<vz_ecoregion2) vz_ecoregion2 = v14ecoregion;

	//---Renumber saf to match numbers in table and shadesets----------
	int vz_ecoregion = con(vz_ecoregion2 == 1, 25, vz_ecoregion2);

	return((PNVcode)vz_ecoregion);

} // end of int PNVmodel::Subregion10010208()


PNVcode PNVmodel::Subregion10010210()
{ //  MBS and GP National Forests, plus MRNP.  Mostly WWC
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		//-----Alpine----------
		float _blalp1 = 140;
		float _blalp2 = 136.5;
		float _blalp3 = -12;
		float _blalp4 = -6;
		float _blalp5 = 135;
		float _blalp6 = 0.5;
		//cad
		float _blalp11 = -.21;
		//fog with PSL not APSL
		float _blalp12 = -100;

		vgtmp2 = int(((( _blalp12 * foggrid + _blalp11 * cadgrid + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + 
				(_blalp5 + _blalp6 * prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -1200;
		float _blpk2 = 150;
		float _blpk3 = -30;
		float _blpk4 = -19;
		float _blpk5 = 150;
		float _blpk6 = 2.0;
		//cad
		float _blpk11 = -.21;
		//fog with PSL not APSL
		float _blpk12 = -100;

		vgtmp3 = int(((( _blpk12 * foggrid + _blpk11 * cadgrid + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.30 * _blpk5 + _blpk6 * prcpgrid) * 
				singrid + (0.70 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		//----Lower MH-----
		float _blmh1 = -1250;
		float _blmh2 = 169;
		float _blmh3 = -315;
		float _blmh4 = -28;
		float _blmh5 = 1.0;
		float _blmh6 = 340;
		float _blmh7 = 3;
		float _blmh8 = -.8;
		float _blmh9 = -.2;
		//cad
		float _blmh11 = -.31;
		//fog with PSL not APSL
		float _blmh12 = -1000;

		vgtmp4 = int((((  _blmh12 * foggrid + _blmh11 * cadgrid + 
							_blmh1 + _blmh2 * msltgrid) + (_blmh3 * pow(tmgrid, 1))) + 
					(_blmh4 + _blmh5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.30 * _blmh6 + _blmh7 * tmgrid) + (_blmh8 + 
					_blmh9 * tmgrid) * prcpgrid) * singrid + 
				(0.70 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -2000;
		float _blsf2 = 148;
		float _blsf3 = -260;
		float _blsf4 = -21.0;
		float _blsf5 = 0.5;
		float _blsf6 = 530;
		float _blsf7 = 10;
		float _blsf8 = -0.5;
		float _blsf9 = -.25;
		//cad
		float _blsf11 = -.50;
		////fog with PSL not APSL
		float _blsf12 = -900;

		vgtmp5 = int(((((_blsf12 * foggrid) + (_blsf11 * cadgrid) + 
							_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.30 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + 
					_blsf9 * tmgrid) * prcpgrid) * singrid + 
				(0.70 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	// create the elevation boundary grid
	{
		float _blwh1 = -17300;
		float _blwh2 = 450;
		float _blwh3 = -350;
		float _blwh4 = -110;
		float _blwh5 = -4.3;
		float _blwh6 = 400;
		float _blwh7 = 30;
		float _blwh8 = 20;
		float _blwh9 = -4;
		//cad
		float _blwh11 = -0.9;
		//fog
		float _blwh12 = -1500;

		vgtmp6 = int(((( _blwh11 * cadgrid + _blwh12 * foggrid + 
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * 
						tmgrid) + (_blwh4 + _blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + 
				 (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 4100;
		float _blsaf2 = 53;
		float _blsaf3 = -225;
		float _blsaf4 = -80;
		float _blsaf5 = 33;
		float _blsaf6 = 220;
		float _blsaf7 = 10;
		float _blsaf8 = -2.5;
		float _blsaf9 = -0.2;
		//cad
		float _blsaf11 = -.31;
		//fog with PSL not APSL
		float _blsaf12 = -200;

		vgtmp7 = int((((_blsaf12 * foggrid + _blsaf11 * cadgrid + 
							_blsaf1 + _blsaf2 * msltgrid) + _blsaf3 * tmgrid) + 
					(_blsaf4 + _blsaf5 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.30 * _blsaf6 + _blsaf7 * tmgrid) + (_blsaf8 + _blsaf9 * 
					tmgrid) * prcpgrid) * singrid + 
				(0.70 * _blsaf6 * (solargrid - 14930) / 3380.4 ));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//*******************************************************************
	{
		//----Upper DF-----
		float _budf1 = 600;
		float _budf2 = -21;
		float _budf3 = -.001;
		float _budf4 = -.00102;
		float _budf5 = 900;
		float _budf6 = -140;
		float _budf7 = 4900;
		float _budf8 = -194;

		vgtmp9 = int(((_budf1 + _budf2 * pow(tmgrid, 2)) + (_budf3 + 
						_budf4 * pow(tmgrid, 2)) * pow(prcpgrid, 3)) + ((_budf5 + _budf6 * 
						pow(tmgrid, .5)) + (_budf7 + _budf8 * tmgrid) / (prcpgrid + 10)) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 700000;
		float _bldf2 = -33900;
		float _bldf3 = 539.1;
		float _bldf4 = -11;
		float _bldf5 = -7;
		float _bldf7 = -1300;
		float _bldf8 = 700;
		float _bldf9 = -10000;
		float _bldf10 = 1900;
		//cad
		float _bldf11 = -0.5;
		//fog with PSL not APSL
		float _bldf12 = 300;
		float _bldf14 = -2.8;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf7 + _bldf8 * pow(tmgrid, .5)) + 
				 (_bldf9 + _bldf10 * tmgrid) / (prcpgrid + 10)) 
				* singrid + (0.60 * _bldf7 * (solargrid - 14930) / 3380.4));

	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		//----Upper MH-----
		float _bumh1 = 8000;
		float _bumh2 = -56;
		float _bumh3 = -69100;
		float _bumh4 = 80;
		float _bumh5 = -900;
		float _bumh6 = 12;
		//cad
		float _bumh11 = 0.2;
		//fog with PSL not APSL
		float _bumh12 = -500;

		vgtmp11 = int((( _bumh12 * foggrid + _bumh11 * cadgrid + _bumh1 + 
						_bumh2 * tmgrid) + ( _bumh3 + 
							_bumh4 * pow(tmgrid, 3)) / prcpgrid) + ((_bumh5 + _bumh6 
							* pow(tmgrid, 2))) * singrid);
	}

	//=====================================================================
	// vgtmp2=alp, vgtmp3=pkl, vgtmp4=mhz, vgtmp5=psfz, vgtmp6=whz, vgtmp7=lbsafz, 
	// vgtmp8=gfz, vgtmp9=ubfz, vgtmp10=lbdfz, vgtmp11=ubmhz,  

	// v10%ecoregion% = con (elevtmp < vgtmp10, 4, 33) /* GL
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 33); // WH
	int v8ecoregion = con (((elevtmp < vgtmp6)), 14, 33); // DF
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp11)), 22, 33); // PSF
	int v6ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp3)), 25, 33); // SAF
	int v4ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

	int vz_ecoregion = v4ecoregion;
	if (v5ecoregion<vz_ecoregion) vz_ecoregion = v5ecoregion;
	if (v6ecoregion<vz_ecoregion) vz_ecoregion = v6ecoregion;
	if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
	if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion;
	if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;

	return((PNVcode)vz_ecoregion);
} // end of Subregion10010210()


PNVcode PNVmodel::Subregion10010103()
{ 
	//
	//----------------------------------------------------------------------
	//            Olympic National Forest/ SW WASH
	//----------------------------------------------------------------------
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	//-----Alpine----------
	{
		float _blalp1 = 240;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;
		//cad
		float _blalp11 = -0.5;
		//fog
		float _blalp12 = -40;

		vgtmp2 = int((((_blalp11 * cadgrid + _blalp12 * foggrid + 
							_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (0.40 * _blalp5 + _blalp6 * prcpgrid) * singrid + 
				(0.60 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PKL) Lower Boundary
	//*******************************************************************
	//-----Parkland-------
	{
		float _blpk1 = -1340;
		float _blpk2 = 150;
		float _blpk3 = -50;
		float _blpk4 = -10;
		float _blpk5 = 150;
		float _blpk6 = 2.0;
		//cad
		float _blpk11 = -0.5;
		//fog
		float _blpk12 = -40;

		vgtmp3 = int((((_blpk11 * cadgrid + _blpk12 * foggrid + 
							_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + 
				(0.40 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.60 * _blpk5 * (solargrid - 14930) / 3380.4));
	}
	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	//----Lower MH-----
	{
		float _blmh1 = -2020;
		float _blmh2 = 186;
		float _blmh3 = -200;
		float _blmh4 = -9;
		float _blmh5 = -3.1;
		float _blmh6 = 450;
		float _blmh7 = -.8;
		float _blmh8 = -.7;
		float _blmh9 = -.2;
		//cad
		float _blmh11 = -.05;
		//fog
		float _blmh12 = -200;

		vgtmp4 = int((((_blmh11 * cadgrid + _blmh12 * foggrid + 
							_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .5)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((_blmh6 + _blmh7 
							* tmgrid) + (0.50 * _blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid + (0.50 * _blmh8 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	//----Lower PSF-----
	{
		float _blsf1 = -3900;
		float _blsf2 = 171;
		float _blsf3 = -290;
		float _blsf4 = -18;
		float _blsf5 = .05;
		float _blsf6 = 230;
		float _blsf7 = 6;
		float _blsf8 = 0.4;
		float _blsf9 = -.35;
		//cad
		//with no CADAELEV
		float _blsf11 = -1.1;
		////fog with PSL not APSL
		float _blsf12 = -450;

		vgtmp5 = int((((_blsf11 * cadgrid + _blsf12 * foggrid + 
							_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + 
				((0.50 * _blsf6 + _blsf7 * tmgrid) + ( _blsf8 + 
					_blsf9 * tmgrid) * prcpgrid) * singrid + 
				(0.50 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	//----Upper SS---new version from ER 104 July 2014--
	{
		float _buss1 = 2200;
		float _buss2 = -115;
		float _buss3 = 40;
		float _buss4 = 250;
		float _buss5 = 1.0;
		float _buss6 = -350;
		float _buss7 = -20;
		float _buss8 = 500;
		float _buss9 = -500;
		//cad
		float _buss12 = 0.7;
		//fog
		float _buss13 = 900;

		vgtmp6 = int((((_buss13 * foggrid + _buss12 * cadgrid + 
							_buss1 + _buss2 * msltgrid) + _buss3 * tmgrid ) + 
					(_buss4 + _buss5 * pow(tmgrid, 0.5)) * 
					pow(prcpgrid, .5)) + ((0.50 * _buss6 + _buss7 * 
						tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) * singrid + 
				(0.50 * _buss6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Western Hemlock (WH) Lower Boundary
	//********************************************************************
	{
		float _blwh1 = -5600;
		float _blwh2 = 150;
		float _blwh3 = -129; 
		float _blwh4 = -40;
		float _blwh5 = -7.5;
		float _blwh6 = 750;
		float _blwh7 = -28;
		float _blwh8 = 60;
		float _blwh9 = -15;
		//CAD
		float _blwh11 = -0.8;
		//FOG
		float _blwh12 = -900.0;

		vgtmp12 = int((((_blwh12 * foggrid + _blwh11 * cadgrid + 
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * tmgrid) + (_blwh4 + 
							_blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//=====================================================================
	//vgtmp2 = lbalp; vgtmp3 = lbpkl; vgtmp4 = lbmhz; vgtmp5 = lbpsfz; vgtmp6=ubssz; vgtmp12 = lbwhz; 

	int v12ecoregion = con (elevtmp < vgtmp12, 16, 33); // GF
	int v10ecoregion = con (elevtmp < vgtmp6, 9, 33); // SS
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 33); // WH
	int v7ecoregion = con (elevtmp < vgtmp4, 22, 33); // SF
	int v6ecoregion = con (elevtmp < vgtmp3, 23, 33); // MH

	int vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	if (v12ecoregion<vz_ecoregion2) vz_ecoregion2 = v12ecoregion;

	return((PNVcode)vz_ecoregion2);

} // end of int PNVmodel::Subregion10010103()


PNVcode PNVmodel::Subregion10010205()
{ // Mt. Hood and N. Willamette NF

	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		float _blalp1 = 240;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;

		vgtmp2 = int((((_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (_blalp5 + _blalp6 * prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -950;
		float _blpk2 = 140; 
		float _blpk3 = -11;
		float _blpk4 = -11;
		float _blpk5 = -200;
		float _blpk6 = 2.1;

		vgtmp3 = int((((_blpk1 + _blpk2 * msltgrid) + _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (_blpk5 + _blpk6 * prcpgrid) * singrid);
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		float _blmh1 = -2700;
		float _blmh2 = 164;
		float _blmh3 = -276;
		float _blmh4 = -6;
		float _blmh5 = -4;
		float _blmh6 = 400;
		float _blmh7 = -1.0;
		float _blmh8 = -.7;
		float _blmh9 = -.2;

		vgtmp4 = int((((_blmh1 + _blmh2 * msltgrid) + (_blmh3 * pow(tmgrid, .5))) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((_blmh6 + _blmh7 
							* tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -1600;
		float _blsf2 = 145;
		float _blsf3 = -370;
		float _blsf4 = -17;
		float _blsf5 = .13;
		float _blsf6 = 480;
		float _blsf7 = -1.8;
		float _blsf8 = -.4;
		float _blsf9 = -.1;

		vgtmp5 = int((((_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + ((_blsf6 + 
							_blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	{
		float _blwh1 = -23200;
		float _blwh2 = 633;
		float _blwh3 = -400;
		float _blwh4 = -74;
		float _blwh5 = -20;
		float _blwh6 = 500;
		float _blwh7 = 100;
		float _blwh8 = 65;
		float _blwh9 = -14;

		vgtmp6 = int((((_blwh1 + _blwh2 * msltgrid) + (_blwh3 * tmgrid)) + (_blwh4 + 
						_blwh5 * tmgrid) * prcpgrid) + ((_blwh6 + _blwh7 * 
							tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 7900;
		float _blsaf2 = -1200;
		float _blsaf3 = -40;
		float _blsaf4 = 10.52;
		float _blsaf5 = 500;
		float _blsaf6 = -27.78;
		float _blsaf7 = -3.19;
		float _blsaf8 = -0.181;

		vgtmp7 = int((( _blsaf1 + _blsaf2 * pow(tmgrid, .5)) + 
					(_blsaf3 + _blsaf4 * pow(tmgrid, .5)) * prcpgrid) + 
				(( _blsaf5 + _blsaf6 * tmgrid) + (_blsaf7 + _blsaf8 * 
					tmgrid) * prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   White Fir (WF) Upper Boundary  
	//*******************************************************************
	{
		float _buwf1 = -12000;
		float _buwf2 = 333;
		float _buwf3 = -270;
		float _buwf4 = -30;
		float _buwf5 = 0.22;
		float _buwf6 = 500;
		float _buwf7 = -5;
		float _buwf8 = -2;
		float _buwf9 = 1.5; 

		vgtmp8 = int((((_buwf1 + _buwf2 * msltgrid) + _buwf3 * tmgrid) + 
					(_buwf4 + _buwf5 * pow(tmgrid, 2)) * prcpgrid) + ((_buwf6 + 
							_buwf7 * tmgrid) + (_buwf8 + _buwf9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   White fir (WF) LOwer Boundary  **LBWF
	//*******************************************************************
	{
		float _blwf1 = -5700;
		float _blwf2 = 200;
		float _blwf3 = -220;
		float _blwf4 = -15;
		float _blwf5 = 0.22;
		float _blwf6 = 800;
		float _blwf7 = -90;
		float _blwf8 = -2;
		float _blwf9 = -0.5;

		vgtmp9 = int((((_blwf1 + _blwf2 * msltgrid) + _blwf3 * tmgrid) + 
					(_blwf4 + _blwf5 * pow(tmgrid, 2)) * prcpgrid) + ((_blwf6 + 
							_blwf7 * tmgrid) + (_blwf8 + _blwf9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 29500;
		float _bldf2 = -538;
		float _bldf3 = 10;
		float _bldf4 = 30;
		float _bldf5 = -16;
		float _bldf6 = -1100;
		float _bldf7 = 200;
		float _bldf8 = -500;
		float _bldf9 = -3000;

		vgtmp10 = int((((_bldf1 + _bldf2 * msltgrid) + _bldf3 * pow(tmgrid, 2)) + 
					(_bldf4 + _bldf5 * tmgrid) * prcpgrid) + (( _bldf6 + _bldf7 * 
							pow(tmgrid, .5)) + (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 10)) 
				* singrid);
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		float _bumh1 = 8000;
		float _bumh2 = -55;
		float _bumh3 = -68000;
		float _bumh4 = 80;
		float _bumh5 = -900;
		float _bumh6 = 12;

		vgtmp11 = int(((_bumh1 + _bumh2 * tmgrid) + ( _bumh3 + 
						_bumh4 * pow(tmgrid, 3)) / prcpgrid) + ((_bumh5 + _bumh6 
						* pow(tmgrid, 2))) * singrid);
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 13;
		float _blgf2 = -1000;
		float _blgf3 = 19200;
		float _blgf4 = 900;
		float _blgf5 = -16;
		float _blgf6 = -21;
		float _blgf7 = 650;
		float _blgf8 = 25;
		float _blgf9 = 65;
		float _blgf10 = -16;

		vgtmp12 = int(((((_blgf1 * pow(msltgrid, 2)) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid) + (_blgf5 + _blgf6 * tmgrid) 
					* prcpgrid) + ((_blgf7 + _blgf8 * tmgrid) + 
						(_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid);
	}

	//=====================================================================
	// vgtmp5 = lbpsf; vgtmp6 = lbwh; vgtmp12 = lbgf; 
	float cadaelev = cadgrid + elevtmp;

	// int v12ecoregion = con (cadaelev < vgtmp10,  5, 33); // steppe
	int v11ecoregion = con (((cadaelev < vgtmp6) and (cadaelev < vgtmp7)), 16, 33); // GF
	int v10ecoregion = con (((cadaelev > vgtmp9) and (cadaelev < vgtmp8)), 20, 33); // WF
	int v9ecoregion = con (cadaelev < vgtmp5, 19, 33); // WH
	int v8ecoregion = con (((cadaelev < vgtmp12) and (cadaelev < vgtmp5) and (cadaelev < vgtmp6)), 14, 33); // DF
	int v7ecoregion = con (((cadaelev < vgtmp4) and (cadaelev < vgtmp11)), 22, 33); // SF
	int v6ecoregion = con (((cadaelev < vgtmp3) and (cadaelev < vgtmp11)), 23, 33); // MH
	int v5ecoregion = con (((cadaelev < vgtmp3) and (cadaelev > vgtmp11) and (cadaelev > vgtmp7)), 25, 33); // SAF
	int v4ecoregion = con (cadaelev < vgtmp2, 32, 33); // PKL

	int vz_ecoregion = v4ecoregion;
	if (v5ecoregion<vz_ecoregion) vz_ecoregion = v5ecoregion;
	if (v6ecoregion<vz_ecoregion) vz_ecoregion = v6ecoregion;
	if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
	if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion;
	if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;
	if (v10ecoregion<vz_ecoregion) vz_ecoregion = v10ecoregion;
	if (v11ecoregion<vz_ecoregion) vz_ecoregion = v11ecoregion;
	// if (v12ecoregion<vz_ecoregion) vz_ecoregion = v12ecoregion;

	return((PNVcode)vz_ecoregion);
} // end of Subregion10010205()


PNVcode PNVmodel::Subregion10010206()
{ // Gifford Pinchot National Forest
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		//-----Alpine----------
		float _blalp1 = 240;
		float _blalp2 = 139;
		float _blalp3 = -12;
		float _blalp4 = -10.028;
		float _blalp5 = 135;
		float _blalp6 = -2.1;
		float _blalp11 = -0.2;
		float _blalp12 = 1;

		vgtmp2 = int((_blalp11 * cadgrid + _blalp12 * foggrid + 
					_blalp1 + _blalp2 * msltgrid + _blalp3 * tmgrid + 
					_blalp4 * prcpgrid) + (0.60 * _blalp5 + _blalp6 * prcpgrid) * singrid + 
				(0.40 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PKL) Lower Boundary
	//*******************************************************************
	{
		//-----Parkland-------
		float _blpk1 = 4000;
		float _blpk2 = 80; 
		float _blpk3 = -100;
		float _blpk4 = -40;
		float _blpk5 = 160;
		float _blpk6 = 2.0;
		float _blpk11 = -0.2;
		float _blpk12 = 100;

		vgtmp3 = int((_blpk11 * cadgrid + _blpk12 * foggrid + 
					_blpk1 + _blpk2 * msltgrid + _blpk3 * tmgrid + 
					_blpk4 * prcpgrid) + (0.40 * _blpk5 + _blpk6 * prcpgrid) * 
				singrid + (0.60 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		float _blmh1 = -900;
		float _blmh2 = 150;
		float _blmh3 = -50;
		float _blmh4 = -30;
		float _blmh5 = -2.0;
		float _blmh6 = 400;
		float _blmh7 = 20;
		float _blmh8 = 10;
		float _blmh9 = -3;
		//cad
		float _blmh11 = -0.5;
		//fog with PSL not APSL
		float _blmh12 = -300;

		vgtmp4 = int((_blmh11 * cadgrid + _blmh12 * foggrid + 
					_blmh1 + _blmh2 * msltgrid + _blmh3 * pow(tmgrid, 1.0) + 
					(_blmh4 + _blmh5 * pow(tmgrid, 1.0)) * prcpgrid) + 
				((0.40 * _blmh6 + _blmh7 * tmgrid) + 
				 (_blmh8 + _blmh9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		float _blsf1 = -1100;
		float _blsf2 = 140;
		float _blsf3 = -40;
		float _blsf4 = -35;
		float _blsf5 = -2.0;
		float _blsf6 = 400;
		float _blsf7 = 40;
		float _blsf8 = 10;
		float _blsf9 = -3;
		//cad
		float _blsf11 = -0.7;
		//fog 
		float _blsf12 = -400;

		vgtmp5 = int((_blsf11 * cadgrid + _blsf12 * foggrid + 
					_blsf1 + _blsf2 * msltgrid + _blsf3 * tmgrid + 
					(_blsf4 + _blsf5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.60 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + 
					_blsf9 * tmgrid) * prcpgrid) * singrid + 
				(0.40 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	{
		float _blwh1 = -11000;
		float _blwh2 = 300;
		float _blwh3 = -150;
		float _blwh4 = -50;
		float _blwh5 = -8;
		float _blwh6 = 550;
		float _blwh7 = 30;
		float _blwh8 = 70;
		float _blwh9 = -14;
		//CAD
		float _blwh11 = -0.6;
		//FOG
		float _blwh12 = -1200.0;

		vgtmp6 = int( _blwh12 * foggrid + _blwh11 * cadgrid + 
				_blwh1 + _blwh2 * msltgrid + _blwh3 * tmgrid + (_blwh4 + _blwh5 * tmgrid) * prcpgrid + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) * 
				singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Lower Boundary
	//*******************************************************************
	{
		//----Lower SAF-------
		float _blsaf1 = 8720;
		float _blsaf2 = -595;
		float _blsaf3 = 20;
		float _blsaf4 = 2;
		float _blsaf5 = 200;
		float _blsaf6 = -27.78;
		float _blsaf7 = -3.19;
		float _blsaf8 = -0.181;
		float _blsaf9 = -60;
		//cad
		float _blsaf11 = -0.7;
		//fog 
		float _blsaf12 = -300;

		vgtmp7 = int(((_blsaf11 * cadgrid + _blsaf12 * foggrid + 
						_blsaf1 + _blsaf9 * msltgrid + _blsaf2 * pow(tmgrid, .5)) + 
					(_blsaf3 + _blsaf4 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.60 * _blsaf5 + _blsaf6 * tmgrid) + (_blsaf7 + _blsaf8 * 
					tmgrid) * prcpgrid) * singrid + 
				(0.40 * _blsaf5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Upper Boundary
	//*******************************************************************
	{
		//----Upper SAF-------
		float _busaf1 = 3800;
		float _busaf2 = 22000;
		float _busaf3 = -180;
		float _busaf4 = -172.6;
		float _busaf5 = 954.5;
		float _busaf6 = -107.2;
		float _busaf7 = -1056.4;
		float _busaf8 = 30000;

		vgtmp8 = int(((_busaf1 + _busaf2 / tmgrid) + 
					(_busaf3 + _busaf4 / tmgrid) * prcpgrid) + 
				(( _busaf5 + _busaf6 * tmgrid) + (_busaf7 + _busaf8 / 
					tmgrid) / prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Upper Boundary
	//*******************************************************************
	{
		float _budf1 = 252700;
		float _budf2 = -9000;
		float _budf3 = 81.7;
		float _budf4 = -120;
		float _budf5 = -90;
		float _budf6 = -9.5;
		float _budf7 = 1000;
		float _budf8 = -80;
		float _budf9 = 20000;
		float _budf10 = -200;
		float _budf11 = -0.7;
		float _budf12 = -300;

		vgtmp9 = int((((_budf11 * cadgrid + _budf12 * foggrid + 
							_budf1 + _budf2 * msltgrid + _budf3 * pow(msltgrid, 2)) + 
						_budf4 * tmgrid) + 
					(_budf5 + _budf6 * tmgrid) * prcpgrid) + ((_budf7 + _budf8 * 
							pow(tmgrid, .5)) + (_budf9 + _budf10 * tmgrid) / (prcpgrid + 10)) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 699800;
		float _bldf2 = -33900;
		float _bldf3 = 539.1;
		float _bldf4 = -11;
		float _bldf5 = -7;
		float _bldf7 = -1300;
		float _bldf8 = 700;
		float _bldf9 = -10000;
		float _bldf10 = 1900;
		//cad
		float _bldf11 = -0.5;
		//fog 
		float _bldf12 = 300;
		float _bldf14 = -2.8;

		vgtmp10 = int( _bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + _bldf14 * pow(msltgrid, 3) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				((0.40 * _bldf7 + _bldf8 * pow(tmgrid, .5)) + 
				 (_bldf9 + _bldf10 * tmgrid) / (prcpgrid + 10)) 
				* singrid + (0.60 * _bldf7 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		//----Upper MH-----
		float _bumh1 = 7700;
		float _bumh2 = -55;
		float _bumh3 = -70000;
		float _bumh4 = 80;
		float _bumh5 = -900;
		float _bumh6 = 12;

		vgtmp11 = int(((_bumh1 + _bumh2 * tmgrid) + ( _bumh3 + 
						_bumh4 * pow(tmgrid, 3)) / prcpgrid) + ((_bumh5 + _bumh6 
						* pow(tmgrid, 2))) * singrid);
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 723300;
		float _blgf2 = -26423;
		float _blgf3 = 241.1;
		float _blgf4 = -300;
		float _blgf5 = 5100;
		float _blgf6 = -150;
		float _blgf7 = 400;
		float _blgf8 = 5;
		float _blgf9 = 25;
		float _blgf10 = -5;
		float _blgf15 = 5;
		//cad
		float _blgf11 = -0.9;
		//fog
		float _blgf12 = -100;

		vgtmp12 = int((_blgf12 * foggrid + _blgf11 * cadgrid + 
					_blgf3 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf1 
					+ _blgf4 * tmgrid + _blgf5 + _blgf6 * prcpgrid + _blgf15 * pow(prcpgrid,2)) + 
				(0.60 * _blgf7 + _blgf8 * tmgrid + 
				 (_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
				(0.40 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary  from 10203  02-18-05
	//*******************************************************************
	{
		float _buss1 = 5700;
		float _buss2 = -250;
		float _buss3 = 200;
		float _buss4 = 160;
		float _buss5 = 25;
		float _buss6 = -10;
		float _buss7 = -0.2;
		float _buss8 = -0.6;
		float _buss9 = -1;
		float _buss10 = -1;

		vgtmp13 = int((((50 * foggrid + 1.41 * cadgrid + _buss1 + _buss2 * msltgrid) + _buss3 * tmgrid + 
						_buss10 * pow(tmgrid, 1)) + (_buss4 + _buss5 * pow(tmgrid, 1)) * pow(prcpgrid, .5)) + 
				((0.60 * _buss6 + _buss7 * 
				  tmgrid) + (_buss8 + _buss9 * tmgrid) * prcpgrid) 
				* singrid + (0.60 * _buss6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*********************************ER 10206**********************************

	//=====================================================================
	//vgtmp2=lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=lbwh(ubgf); vgtmp7=lbsaf;
	//vgtmp8=ubsaf; vgtmp9=ubdf; vgtmp10=lbdf; vgtmp11=ubmh; vgtmp12=lbgf; vgtmp13=ubss

	// int v12ecoregion = con (((elevtmp < vgtmp6) and (elevtmp < vgtmp10) and (elevtmp < vgtmp12)), 4, 33); // GL, 
	int v11ecoregion = con (elevtmp < vgtmp13, 9, 33); // SS
	int v10ecoregion = con (((elevtmp > vgtmp12) and (elevtmp < vgtmp6)), 16, 33); // GF
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 33); // WH
	int v8ecoregion = con (((elevtmp < vgtmp6) and (elevtmp < vgtmp12) /* and (elevtmp > vgtmp10) */), 14, 33); // DF
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp11)), 22, 33); // SF
	int v6ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp11) and (elevtmp > vgtmp7)), 25, 33); // SAF
	int v4ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

	int vz_ecoregion = v4ecoregion;
	if (v5ecoregion<vz_ecoregion) vz_ecoregion = v5ecoregion;
	if (v6ecoregion<vz_ecoregion) vz_ecoregion = v6ecoregion;
	if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
	if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion;
	if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;
	if (v10ecoregion<vz_ecoregion) vz_ecoregion = v10ecoregion;
	if (v11ecoregion<vz_ecoregion) vz_ecoregion = v11ecoregion;

	return((PNVcode)vz_ecoregion);
} // end of Subregion10010206()


PNVcode PNVmodel::Subregion10010214()
{ // N. Willamette Valley, with a little bit in WWC
	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		//----Lower MH-----
		float _blmh1 = -4500;
		float _blmh2 = 209;
		float _blmh3 = -234;
		float _blmh4 = -13;
		float _blmh5 = -4;
		float _blmh6 = 450;
		float _blmh7 = 2;
		float _blmh8 = -.5;
		float _blmh9 = -.2;
		//cad
		float _blmh11 = -0.1;
		//fog
		float _blmh12 = -1000;

		vgtmp4 = int((((_blmh11 * cadgrid + _blmh12 * foggrid + 
							_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .5)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((_blmh6 + _blmh7 
							* tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (PSF) Lower Boundary
	//********************************************************************
	{
		//----Lower PSF-----
		float _blsf1 = -5550;
		float _blsf2 = 251;
		float _blsf3 = -430;
		float _blsf4 = -26;
		float _blsf5 = .22;
		float _blsf6 = 400;
		float _blsf7 = 6;
		float _blsf8 = .4;
		float _blsf9 = -.35;
		//cad
		float _blsf11 = -0.1;
		//fog
		float _blsf12 = -2000;

		vgtmp5 = int((((_blsf11 * cadgrid + _blsf12 * foggrid + 
							_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + ((_blsf6 + 
							_blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	{
		//----Upper SS-----
		float _buss1 = 14300;
		float _buss2 = -445;
		float _buss3 = 120;
		float _buss4 = 625;
		float _buss5 = 40;
		float _buss6 = -300;
		float _buss7 = -20;
		float _buss8 = -2000;
		float _buss9 = -500;

		vgtmp6 = int((((_buss1 + _buss2 * msltgrid) + _buss3 * tmgrid) + 
					(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5)) + ((_buss6 + _buss7 * 
						tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Western Hemlock Lower Boundary (UBDF) 
	//*******************************************************************
	{
		//----LBWH / Upper DF...  -----
		float _blwh1 = -29800;
		float _blwh2 = 690;
		float _blwh3 = -239;
		float _blwh4 = -60;
		float _blwh5 = -15.5;
		float _blwh6 = 150;
		float _blwh7 = 75;
		float _blwh8 = 55;
		float _blwh9 = -10;
		//cad
		float _blwh11 = -0.5;
		//fog
		float _blwh12 = -1100;

		vgtmp9 = int(((_blwh11 * cadgrid + _blwh12 * foggrid + 
						_blwh1 + _blwh2 * msltgrid + _blwh3 * tmgrid) + (_blwh4 + 
							_blwh5 * tmgrid) * prcpgrid) + ((_blwh6 + _blwh7 * 
								tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		//----Lower DF-----
		float _bldf1 = -7600;
		float _bldf2 = 220;
		float _bldf3 = -60;
		float _bldf4 = -54;
		float _bldf5 = -4;
		float _bldf6 = 600;
		float _bldf7 = 200;
		float _bldf8 = -500;
		float _bldf9 = -3000;
		//cad
		float _bldf11 = -0.1;
		//fog
		float _bldf12 = -2500;

		vgtmp10 = int((((_bldf11 * cadgrid + _bldf12 * foggrid + 
							_bldf1 + _bldf2 * msltgrid) + _bldf3 * pow(tmgrid, 2)) + 
					(_bldf4 + _bldf5 * tmgrid) * prcpgrid) + (( _bldf6 + _bldf7 * 
							pow(tmgrid, .5)) + (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 10)) * singrid);
	}

	//********************************************************************
	//    GRAND FIR (GF) Upper Boundary
	//********************************************************************
	{
		//----Upper GF-----new for ER9e
		//float _bugf1 -2200
		float _bugf1 = -2350;
		// float _bugf2 = 175;
		float _bugf3 = 16;
		float _bugf4 = -.003;
		float _bugf5 = -.00058;
		float _bugf6 = 1000;
		float _bugf7 = -140;
		float _bugf8 = 5000;
		float _bugf9 = -194;
		float _bugf11 = -88800;
		float _bugf12 = 3205;
		float _bugf13 = -27.95;

		vgtmp12 = int((((_bugf1 + (_bugf11 + _bugf12 * msltgrid + _bugf13 * pow(msltgrid, 2))) + 
						_bugf3 * pow(tmgrid, 2)) + (_bugf4 + 
						_bugf5 * pow(tmgrid, 2)) * pow(prcpgrid, 2)) + ((_bugf6 + _bugf7 * 
						pow(tmgrid, .5)) + (_bugf8 + _bugf9 * tmgrid) / (prcpgrid + 10)) 
				* singrid);
	} 

	//********************************************************************
	//    GRAND FIR (GF) Lower Boundary
	//********************************************************************
	{
		//----Lower GF-----new for ER9e
		float _blgf1 = 7450;
		float _blgf2 = 18;
		float _blgf3 = -150;
		float _blgf4 = -88;
		float _blgf5 = -10;
		float _blgf6 = 500;
		float _blgf7 = -90;
		float _blgf8 = 28;
		float _blgf9 = -1;
		////cad
		float _blgf11 = -0.6;
		////fog
		float _blgf12 = -1000;

		vgtmp13 = int(((_blgf11 * cadgrid + _blgf12 * foggrid + 
						_blgf1 + _blgf2 * msltgrid) + 
					(_blgf3 * tmgrid) + (_blgf4 + _blgf5 * tmgrid) * prcpgrid) + 
				((_blgf6 + _blgf7 * tmgrid) + (_blgf8 + _blgf9 * tmgrid) * prcpgrid) * singrid);
	}

	//======================================== 10214 ============================
	//vgtmp2=lbalp, 3=lbpkl 4=lbmh, 5=lbpsf, 6=ubss, 
	//9=lbwh, 10=lbdf, 11=ubmh, 12=ubgf,  13=lbgf

	int v12ecoregion = con (((elevtmp > vgtmp13) and (elevtmp > vgtmp10) and (elevtmp < vgtmp12)), 16, 22); // GF
	int v10ecoregion = con (elevtmp < vgtmp6, 9, 22); // SS
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 22); // WH
	int v8ecoregion = con (elevtmp < vgtmp9, 18, 22); // DF
	int v7ecoregion = con ((elevtmp < vgtmp4), 22, 22); // PSF
	int v6ecoregion = con (elevtmp < vgtmp10, 4, 22); // W Grassland

	int vz_ecoregion = v6ecoregion;
	if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
	if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion;
	if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;
	if (v10ecoregion<vz_ecoregion) vz_ecoregion = v10ecoregion;
	if (v12ecoregion<vz_ecoregion) vz_ecoregion = v12ecoregion;

	return((PNVcode)vz_ecoregion);

} // end of Subregion10010214()


PNVcode PNVmodel::Subregion10010215()
{ 
	//----------------------------------------------------------------------
	//            Puget Lowlands,  
	//----------------------------------------------------------------------

	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALP) Lower Boundary
	//*******************************************************************
	{
		//-----Alpine----------
		float _blalp1 = -250;
		float _blalp2 = 137;
		float _blalp3 = -15;
		float _blalp4 = -5;
		float _blalp5 = 135;
		float _blalp6 = 0.5;

		vgtmp2 = int((((_blalp1 + _blalp2 * msltgrid) + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + (_blalp5 + _blalp6 * prcpgrid) * singrid);

	}
	//*******************************************************************
	//*******************************************************************
	//   Parkland (PK) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -1400;
		float _blpk2 = 150;
		float _blpk3 = -50;
		float _blpk4 = -10;
		float _blpk5 = 180;
		float _blpk6 = 2;

		vgtmp3 = int((((_blpk1 + _blpk2 * msltgrid) +  _blpk3 * tmgrid) + 
					(_blpk4 * prcpgrid)) + (0.50 + _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Lower Boundary
	//********************************************************************
	{
		//----Lower MH-----
		float _blmh1 = -6350;
		float _blmh2 = 252;
		float _blmh3 = -500;
		float _blmh4 = -12;
		float _blmh5 = -2.5;
		float _blmh6 = 360;
		float _blmh7 = -.8;
		float _blmh8 = -.7;
		float _blmh9 = -.2;

		vgtmp4 = int((((_blmh1 + _blmh2 * msltgrid) + _blmh3 * pow(tmgrid, .6)) + 
					(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + ((_blmh6 + _blmh7 
							* tmgrid) + (_blmh8 + _blmh9 * tmgrid) * prcpgrid) 
				* singrid);
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (SF) Lower Boundary
	//********************************************************************
	{
		//----Lower PSF-----
		float _blsf1 = -5400;
		float _blsf2 = 230;
		float _blsf3 = -390;
		float _blsf4 = -28;
		float _blsf5 = .15;
		float _blsf6 = 700;
		float _blsf7 = 3;
		float _blsf8 = -.5;
		float _blsf9 = -.10;
		//CAD
		float _blsf11 = -0.2;
		//FOG
		float _blsf12 = -400.0;

		vgtmp5 = int((((_blsf12 * foggrid + _blsf11 * cadgrid + 
							_blsf1 + _blsf2 * msltgrid) + _blsf3 * tmgrid) + 
					(_blsf4 + _blsf5 * pow(tmgrid, 2)) * prcpgrid) + ((0.50 * _blsf6 + 
							_blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * prcpgrid) 
				* singrid + (0.50 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Sitka Spruce (SS) Upper Boundary
	//*******************************************************************
	{
		//----Upper SS-----
		float _buss1 = 7000;
		float _buss2 = -264;
		float _buss3 = 77;
		float _buss4 = 200;
		float _buss5 = 80;
		float _buss6 = -80;
		float _buss7 = -20;
		float _buss8 = 500;
		float _buss9 = -500;

		vgtmp6 = int((((_buss1 + _buss2 * msltgrid) + _buss3 * tmgrid) + 
					(_buss4 + _buss5 * tmgrid) * pow(prcpgrid, .5)) + ((_buss6 + _buss7 * 
						tmgrid) + (_buss8 + _buss9 / tmgrid) / prcpgrid) 
				* singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	{
		float _blwh1 = -7300;
		float _blwh2 = 202;
		float _blwh3 = -100;
		float _blwh4 = -30;
		float _blwh5 = -8;
		float _blwh6 = 600;
		float _blwh7 = 40;
		float _blwh8 = 70;
		float _blwh9 = -18;
		//CAD
		float _blwh11 = -0.5;
		//FOG
		float _blwh12 = -600.0;

		vgtmp7 = int(((( _blwh12 * foggrid + _blwh11 * cadgrid + 
							_blwh1 + _blwh2 * msltgrid) + _blwh3 * tmgrid) + (_blwh4 + 
							_blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + 
					_blwh9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAF) Upper Boundary
	//*******************************************************************
	{
		//----Upper SAF-------
		float _busaf1 = 4000;
		float _busaf2 = 22000;
		float _busaf3 = -180;
		float _busaf4 = -172.6;
		float _busaf5 = 954.5;
		float _busaf6 = -107.2;
		float _busaf7 = -1056.4;
		float _busaf8 = 30000;

		vgtmp8 = int(((_busaf1 + _busaf2 / tmgrid) + 
					(_busaf3 + _busaf4 / tmgrid) * prcpgrid) + 
				(( _busaf5 + _busaf6 * tmgrid) + (_busaf7 + _busaf8 / 
					tmgrid) / prcpgrid) * singrid);
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 232600;
		float _bldf2 = -8825;
		float _bldf3 = 81.57;
		float _bldf4 = 30;
		float _bldf5 = 4;
		float _bldf6 = 5;
		float _bldf7 = -100;
		float _bldf8 = -10;
		float _bldf9 = -5000;
		//CAD
		float _bldf11 = -0.4;
		//FOG
		float _bldf12 = 400.0;

		vgtmp10 = int(_bldf12 * foggrid + _bldf11 * cadgrid + 
				_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + 
				(_bldf4 + _bldf5 * tmgrid) * prcpgrid + 
				(0.40 * _bldf6 + _bldf7 * pow(tmgrid, .5) + (_bldf8 + _bldf9 * tmgrid) / (prcpgrid + 0)) 
				* singrid + (0.60 * _bldf6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		//----Upper MH-----
		float _bumh1 = 7500;
		float _bumh2 = -90;
		float _bumh3 = -69000;
		float _bumh4 = 80;
		float _bumh5 = -900;
		float _bumh6 = 12;

		vgtmp11 = int(((_bumh1 + _bumh2 * tmgrid) + ( _bumh3 + 
						_bumh4 * pow(tmgrid, 3)) / prcpgrid) + ((_bumh5 + _bumh6 
						* pow(tmgrid, 2))) * singrid);
	}

	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf1 = 364.12;
		float _blgf2 = -40175.5;
		float _blgf3 = 1113300;
		float _blgf4 = -700;
		float _blgf5 = -125;
		float _blgf6 = 2;
		float _blgf7 = 450;
		float _blgf8 = 10;
		float _blgf9 = 25;
		float _blgf10 = -6;
		//cad
		float _blgf11 = -0.4;
		//fog
		float _blgf12 = -500;
		float _blgf14 = -100;
		float _blgf15 = 3;

		vgtmp12 = int(((((_blgf12 * foggrid) + (_blgf11 * cadgrid) + 
							_blgf1 * pow(msltgrid, 2) + _blgf2 * msltgrid + _blgf3) 
						+ _blgf4 * tmgrid) + (_blgf5 + _blgf6 * pow(tmgrid, 2)) 
					* prcpgrid) + _blgf14 + _blgf15 * pow(msltgrid, 2) + 
				((0.5 * _blgf7 + _blgf8 * tmgrid) + 
				 (_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
				(0.50 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//=====================================================================
	//vgtmp2=lbalp; vgtmp3=lbpkl; vgtmp4=lbmh; vgtmp5=lbpsf; vgtmp6=ssz; vgtmp7=lbwh;
	//vgtmp8=ubsaf; vgtmp9=ubdf; vgtmp10=lbdf; vgtmp11=ubmh; vgtmp12=lbgf;  

	// v12%ecoregion% = con (( (elevtmp < vgtmp12) and (elevtmp < vgtmp10) and (elevtmp < vgtmp7)), 4, 33) /* W Grassl
	int v11ecoregion = con (((elevtmp < vgtmp7) and (elevtmp < vgtmp5) ), 16, 33); // GF
	int v10ecoregion = con (elevtmp < vgtmp6, 9, 33); // SS
	int v9ecoregion = con (elevtmp < vgtmp5, 19, 33); // WH
	int v8ecoregion = con (( (elevtmp < vgtmp7) and (elevtmp < vgtmp12) ), 14, 33); // DF
	int v7ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp11)), 22, 33); // SF
	int v6ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
	int v5ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp11) and (elevtmp > vgtmp7)), 25, 33); // SAF
	int v4ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

	// vz_ecoregion = min( v4ecoregion, v5ecoregion, v6ecoregion, v7ecoregion, 
	//     v8ecoregion, v9ecoregion, v10ecoregion, v11ecoregion   ) 

	int vz_ecoregion2 = v4ecoregion;
	if (v5ecoregion<vz_ecoregion2) vz_ecoregion2 = v5ecoregion;
	if (v6ecoregion<vz_ecoregion2) vz_ecoregion2 = v6ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v8ecoregion<vz_ecoregion2) vz_ecoregion2 = v8ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;
	// if (v12ecoregion<vz_ecoregion2) vz_ecoregion2 = v12ecoregion;

	return((PNVcode)vz_ecoregion2);
} // end of int PNVmodel::Subregion10010215()


PNVcode PNVmodel::Subregion10020101()
{ // Okanogan National Forest

	//*************************ecoregion 10020101******************************************
	//*******************************************************************
	//   Alpine (ALPZ) Lower Boundary
	//*******************************************************************
	{
		float _blalp1 = -380;
		float _blalp2 = 146;
		float _blalp3 = -15;
		float _blalp4 = -9;
		float _blalp5 = 105; 
		float _blalp6 = -1.5;

		vgtmp2 = int(((_blalp1 + _blalp2 * msltgrid + _blalp3 * tmgrid) + 
					(_blalp4 * prcpgrid)) + 
				(0.50 * _blalp5 + _blalp6 * prcpgrid) * singrid + 
				(0.50 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PKZ) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -1400;
		float _blpk2 = 149;
		float _blpk3 = -22;
		float _blpk4 = -7;
		float _blpk5 = 250;
		float _blpk6 = -1.5;
		float _blpk11 = -0.3;

		vgtmp3 = int((((_blpk11 * cadgrid + _blpk1 + _blpk2 * msltgrid) + 
						(_blpk3 * tmgrid)) + (_blpk4 * prcpgrid)) + 
				(0.60 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.40 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAFZ) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 4450;
		float _blsaf2 = 95;
		float _blsaf3 = -400;
		float _blsaf4 = -30;
		float _blsaf5 = 9.0;
		float _blsaf6 = 600;
		float _blsaf7 = 10;
		float _blsaf8 = -2.5;
		float _blsaf9 = -1;
		float _blsaf11 = -0.6; 

		vgtmp4 = int((((_blsaf11 * cadgrid + _blsaf1 + _blsaf2 * msltgrid) + 
						(_blsaf3 * pow(tmgrid, 1))) + 
					(_blsaf4 + _blsaf5 * pow(tmgrid, .5)) * prcpgrid) + 
				((0.50 * _blsaf6 + _blsaf7 * tmgrid) + (_blsaf8 + _blsaf9 * 
					tmgrid) * prcpgrid) * singrid + (0.50 * _blsaf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _lbdf1 = 253500;
		float _lbdf2 = -9000;
		float _lbdf3 = 80.7;
		float _lbdf4 = -100;
		float _lbdf5 = -70;
		float _lbdf6 = -4.4;
		float _lbdf7 = 500;
		float _lbdf8 = -55;
		float _lbdf9 = 5000;
		float _lbdf10 = -195;

		float _bldf11 = -0.5;
		float _bldf12 = -1000;

		vgtmp5 = int((((_bldf11 * cadgrid  + _bldf12 * foggrid + _lbdf1 + _lbdf2 * msltgrid) + 
						_lbdf3 * pow(msltgrid, 2) + (_lbdf4 * tmgrid)) + 
					(_lbdf5 + _lbdf6 * tmgrid) * prcpgrid) + ((0.50 * _lbdf7 + _lbdf8 * 
							pow(tmgrid, .5)) + (_lbdf9 + _lbdf10 * tmgrid) / prcpgrid) 
				* singrid + (0.50 * _lbdf7 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Lower PSFZ (SF) / Western Hemlock (WH) Upper Boundary 
	//*******************************************************************
	{
		float _blsf1 = -4000;
		float _blsf2 = 260;
		float _blsf3 = -400;
		float _blsf4 = -50;
		float _blsf5 = -4.0;
		float _blsf6 = 500;
		float _blsf7 = 40;
		float _blsf8 = 15;
		float _blsf9 = -3;
		//cad
		float _blsf11 = -1.0;
		//fog
		float _blsf12 = -1200;

		vgtmp6 = int((((_blsf11 * cadgrid  + _blsf12 * foggrid + 
							_blsf1 + _blsf2 * msltgrid) + (_blsf3 * tmgrid)) + 
					(_blsf4 + _blsf5 * tmgrid) * prcpgrid) + 
				((0.40 * _blsf6 + _blsf7 * tmgrid) + (_blsf8 + _blsf9 * tmgrid) * 
				 prcpgrid) * singrid + (0.60 * _blsf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock Lower Boundary
	//*******************************************************************
	{
		float _blwh1 = -600;
		float _blwh2 = 195;
		float _blwh3 = -420;  
		float _blwh4 = -70;
		float _blwh5 = -4.8;
		float _blwh6 = 450;
		float _blwh7 = 40;
		float _blwh8 = 95;
		float _blwh9 = -18;
		float _blwh11 = -1.5;
		float _blwh12 = -1000;

		vgtmp7 = int((((_blwh11 * cadgrid + _blwh12 * foggrid + _blwh1 + _blwh2 * msltgrid) + 
						(_blwh3 * tmgrid)) + (_blwh4 + _blwh5 * tmgrid) * prcpgrid) + 
				((0.50 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + _blwh9 * tmgrid) * 
				 prcpgrid) * singrid + (0.50 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		float _bumh1 = 19700;
		float _bumh2 = -260;
		float _bumh3 = 170;
		float _bumh4 = -80000;
		float _bumh5 = 19;
		float _bumh6 = -1300;
		float _bumh7 = -15;
		float _bumh11 = 2.0;

		vgtmp8 = int(((( _bumh11 * cadgrid + _bumh1 + _bumh2 * msltgrid) + 
						( _bumh3 * tmgrid)) + ( _bumh4 + _bumh5 * pow(tmgrid, 3)) / prcpgrid) + 
				(0.40 * _bumh6 + _bumh7 * pow(tmgrid, 2)) * singrid + 
				(0.60 * _bumh6 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    LBMHZ / old Pacific Silver Fir (PSF) Upper Boundary 
	//********************************************************************
	{
		float _blmh1 = -1350;
		float _blmh2 = 215;
		float _blmh3 = -350;
		float _blmh4 = -39;
		float _blmh5 = -2.5;
		float _blmh6 = 650; 
		float _blmh7 = -30;     
		float _blmh8 = 1.5;
		float _blmh9 = -.1;
		float _blmh11 = -1.0;
		float _blmh12 = -1400;

		vgtmp9 = int((((_blmh11 * cadgrid  + _blmh12 * foggrid + 
							_blmh1 + _blmh2 * msltgrid) + (_blmh3 * pow(tmgrid, 1))) +  
					(_blmh4 + _blmh5 * pow(tmgrid, 1)) * prcpgrid) + 
				((0.40 * _blmh6 + _blmh7 * tmgrid) + (_blmh8 + _blmh9 * tmgrid) * 
				 prcpgrid) * singrid + 
				(0.60 * _blmh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//********************************************************************
	//    Ponderosa Pine (PP) Lower Boundary
	//********************************************************************
	{
		float _blpp1 = 228900;
		float _blpp2 = -7980;
		float _blpp3 = 69.9;
		float _blpp4 = -120;
		float _blpp5 = -77;
		float _blpp6 = -4.3;
		float _blpp7 = 300;
		float _blpp8 = 20;
		float _blpp9 = 80;
		float _blpp10 = -17;
		float _blpp11 = -0.8;
		float _blpp12 = -1000;

		vgtmp11 = int((((_blpp11 * cadgrid  + _blpp12 * foggrid + _blpp1 + _blpp2 * msltgrid) + 
						_blpp3 * pow(msltgrid, 2) + (_blpp4 * tmgrid   )) + ((_blpp5 + 
							_blpp6 * tmgrid) * prcpgrid) + (0.40 * _blpp7 * _blpp8 * tmgrid) + 
					( _blpp9 + _blpp10 * tmgrid)  / prcpgrid) * singrid + 
				(0.60 * _blpp6 * (solargrid - 14930) / 3380.4));
	}

	//******************************************************************** 
	//vgtmp 2=lbalp, vgtmp3=lbpkl, vgtmp4=lbsafz, vgtmp5=lbdfz, vgtmp6=ubwhz(lbpsfz), 
	//vgtmp7=lbwhz, vgtmp8=ubmhz, vgtmp9=ubpsfz/lbmhz, vgtmp10=lbgfz, vgtmp11=lbppz

	int v11ecoregion = con (vgtmp5 > elevtmp, 10, 33); // PPZ
	// int v9ecoregion = con (((vgtmp5 > elevtmp) and (vgtmp11 > elevtmp)), 5, 33); // ST  
	int v7ecoregion = con (((vgtmp7 > elevtmp)  and (vgtmp3 > elevtmp) and (vgtmp4 > elevtmp)), 14, 33); // DF 
	int v6ecoregion = con (((vgtmp6 > elevtmp) and (vgtmp8 > elevtmp)), 19, 33); // WH
	int v5ecoregion = con (((vgtmp9 > elevtmp) and (vgtmp6 < elevtmp) and (vgtmp3 > elevtmp) and (vgtmp8 > elevtmp)), 22, 33); // PSF 
	int v4ecoregion = con (((vgtmp3 > elevtmp) and (vgtmp9 < elevtmp) and (vgtmp8 > elevtmp) ), 23, 33); // MH 
	int v3ecoregion = con (((vgtmp3 > elevtmp) and (vgtmp8 < elevtmp)), 25, 33); // SAF  
	int v2ecoregion = con (vgtmp2 > elevtmp, 32, 33); //  Pkl

	// Note:  Ponderosa Pine is included  
	int vz_ecoregion = v2ecoregion;
	if(v3ecoregion<vz_ecoregion) vz_ecoregion = v3ecoregion;
	if(v4ecoregion<vz_ecoregion) vz_ecoregion = v4ecoregion;
	if(v5ecoregion<vz_ecoregion) vz_ecoregion = v5ecoregion;
	if(v6ecoregion<vz_ecoregion) vz_ecoregion = v6ecoregion;
	if(v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
	// if(v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;
	if(v11ecoregion<vz_ecoregion) vz_ecoregion = v11ecoregion;

	return((PNVcode)vz_ecoregion);
} // end of Subregion1020101()

/* 10/30/14 use 10020101  
   PNVcode PNVmodel::Subregion10020110()
   { // Far north Wenatchee N.F./ Okanogan N.F.

// *************************ecoregion 10020110 ******************************************
// *******************************************************************
//   Alpine (ALPZ) Lower Boundary
// *******************************************************************
{
float _blalp1 = 50;
float _blalp2 = 145;
float _blalp3 = -12;
float _blalp4 = -8;
float _blalp5 = 250;
float _blalp6 = -2;

vgtmp2 = int((((_blalp1 + _blalp2 * msltgrid) + (_blalp3 * tmgrid)) + 
(_blalp4 * prcpgrid)) + 
(0.50 * _blalp5 + _blalp6 * prcpgrid) * singrid + (0.50 * _blalp5 * (solargrid - 14930) / 3380.4));
}
// *******************************************************************
// *******************************************************************
//   Parkland (PKZ) Lower Boundary
// *******************************************************************
{
float _blpk1 = -1110;
float _blpk2 = 149;
float _blpk3 = -52;
float _blpk4 = -10;
float _blpk5 = 130;
float _blpk6 = 2;
float _blpk11 = -0.3;
//fog
float _blpk12 = -350;

vgtmp3 = int((((  _blpk11 * cadgrid + _blpk12 * foggrid + _blpk1 + _blpk2 * msltgrid) + (_blpk3 * tmgrid)) + 
(_blpk4 * prcpgrid)) + 
(0.50 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
(0.50 * _blpk5 * (solargrid - 14930) / 3380.4));
}

// *******************************************************************
// *******************************************************************
//   Mountain Hemlock (MH) Lower Boundary
// *******************************************************************
{
float _blmh1 = -1150;
float _blmh2 = 197;
float _blmh3 = -330;
float _blmh4 = -30;
float _blmh5 = -4.9;
float _blmh6 = 550;
float _blmh7 = -30;
float _blmh8 = 1.5;
float _blmh9 = -0.1;
float _blmh11 = -2.4;
//fog
float _blmh12 = -650;

vgtmp4 = int((((_blmh11 * cadgrid + _blmh12 * foggrid +  _blmh1 + _blmh2 * 
msltgrid) + (_blmh3 * pow(tmgrid, .5))) + 
(_blmh4 + _blmh5 * pow(tmgrid, .5)) * prcpgrid) + 
(( _blmh6 + _blmh7 * tmgrid) + (0.40 * _blmh8 + _blmh9 * 
tmgrid) * prcpgrid) * singrid + (0.65 * _blmh8 * (solargrid - 14930) / 3380.4));
}

// ********************************************************************
// ********************************************************************
//    Pacific Silver Fir (PSF) Lower Boundary
// ********************************************************************
{
float _blsf1 = -4100;
float _blsf2 = 236;
float _blsf3 = -300;
float _blsf4 = -39;
float _blsf5 = -4.7;
float _blsf6 = 650;
float _blsf7 = 50;
float _blsf8 = -12;
float _blsf9 = 2;
float _blsf11 = -1.2;
//fog
float _blsf12 = -1000;

vgtmp5 = int((((_blsf11 * cadgrid + _blsf12 * foggrid + _blsf1 + _blsf2 * 
					msltgrid) + (_blsf3 * tmgrid)) + 
			(_blsf4 + _blsf5 * tmgrid) * prcpgrid) + 
		((_blsf6 + _blsf7 * tmgrid) + (0.40 * _blsf8 + _blsf9 * tmgrid) * 
		 prcpgrid) * singrid + (0.60 * _blsf8 * (solargrid - 14930) / 3380.4));
}

// *******************************************************************
// *******************************************************************
//    Western Hemlock (WH) Lower Boundary 
// *******************************************************************
{
	float _blwh1 = -6580;
	float _blwh2 = 305;
	float _blwh3 = -540;
	float _blwh4 = -31;
	float _blwh5 = -7;
	float _blwh6 = 1200;
	float _blwh7 = 115;
	float _blwh8 = 60;
	float _blwh9 = -14;
	//cad
	float _blwh11 = -1.2;
	//fog
	float _blwh12 = -1000;

	vgtmp6 = int(((( _blwh11 * cadgrid + _blwh12 * foggrid + _blwh1 + _blwh2 * 
						msltgrid) + (_blwh3 * tmgrid)) + (_blwh4 + 
						_blwh5 * tmgrid) * prcpgrid) + ((0.50 * _blwh6 + _blwh7 * 
							tmgrid) + (_blwh8 + _blwh9 * tmgrid) * prcpgrid) 
			* singrid + (0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
}

// *******************************************************************
// *******************************************************************
//   Subalpine Fir (SAFZ) Lower Boundary
// *******************************************************************
{ 
	float _blsaf1 = 2400;
	float _blsaf2 = 95;
	float _blsaf3 = -700;
	float _blsaf4 = -20;
	float _blsaf5 = -10;
	float _blsaf6 = 500;
	float _blsaf7 = 25;
	float _blsaf8 = -2.5;
	float _blsaf9 = -0.2;
	float _blsaf11 = -0.9; 
	//fog
	float _blsaf12 = -400;

	vgtmp7 = int((((  _blsaf11 * cadgrid + _blsaf12 * foggrid +  _blsaf1 + _blsaf2 * 
						msltgrid) + (_blsaf3 * pow(tmgrid, .5))) + 
				(_blsaf4 + _blsaf5 * pow(tmgrid, .5)) * prcpgrid) + 
			(( 0.40 * _blsaf6 + _blsaf7 * tmgrid) + (_blsaf8 + _blsaf9 * 
				tmgrid) * prcpgrid) * singrid + 
			(0.60 * _blsaf6 * (solargrid - 14930) / 3380.4));
}

// *******************************************************************
// *******************************************************************
//   Douglas-fir (DF) Lower Boundary
// *******************************************************************
{
	float _bldf1 = 254100;
	float _bldf2 = -9000;
	float _bldf3 = 80.7;
	float _bldf4 = -100;
	float _bldf5 = -83;
	float _bldf6 = -4.4;
	float _bldf7 = 480;
	float _bldf8 = -55;
	float _bldf9 = 5000;
	float _bldf10 = -195;
	float _bldf11 = -0.5;
	//fog
	float _bldf12 = -1000;

	vgtmp9 = int((((  _bldf11 * cadgrid + _bldf12 * foggrid + _bldf1 + _bldf2 * 
						msltgrid + _bldf3 * pow(msltgrid, 2)) + (_bldf4 * 
						tmgrid)) + (_bldf5 + _bldf6 * tmgrid) * prcpgrid) + 
			((0.50 * _bldf7 + _bldf8 * pow(tmgrid, .5)) + ( _bldf9 + 
				_bldf10 * tmgrid) / prcpgrid) * singrid + 
			(0.50 * _bldf7 * (solargrid - 14930) / 3380.4));
}

// ********************************************************************
// ********************************************************************
//    Ponderosa Pine (PP) Lower Boundary
// ********************************************************************
{
	float _blpp1 = 230400;
	float _blpp2 = -7980;
	float _blpp3 = 69.9;
	float _blpp4 = -120;
	float _blpp5 = -84;
	float _blpp6 = -4.3;
	float _blpp7 = 300;
	float _blpp8 = 20;
	float _blpp9 = 80;
	float _blpp10 = -17;
	float _blpp11 = -0.8;
	//fog
	float _blpp12 = -1000;

	vgtmp10 = int(((( _blpp11 * cadgrid + _blpp12 * foggrid + _blpp1 + _blpp2 * 
						msltgrid + _blpp3 * pow(msltgrid, 2) ) + (_blpp4 * 
						tmgrid)) + (_blpp5 + _blpp6 * tmgrid) * prcpgrid) + 
			((0.50 * _blpp7 + _blpp8 * tmgrid) + ( _blpp9 + _blpp10 
				* tmgrid) / prcpgrid) * singrid + 
			(0.50 * _blpp7 * (solargrid - 14930) / 3380.4));
}

// ********************************************************************
// ********************************************************************
//    Mountain Hemlock (MH) Upper Boundary
// ********************************************************************
{
	float _bumh1 = 12940;
	float _bumh2 = -118;
	float _bumh3 = -26;
	float _bumh4 = -69500;
	float _bumh5 = 755;
	float _bumh6 = -700;
	float _bumh7 = 14;
	// cad
	float _bumh11 = 1.2;
	//fog
	float _bumh12 = 200;

	vgtmp11 = int(((( _bumh11 * cadgrid + _bumh12 * foggrid + _bumh1 + _bumh2 * msltgrid) + ( _bumh3 * tmgrid)) + ( _bumh4 + 
					_bumh5 * pow(tmgrid, 2)) / prcpgrid) + (0.40 * _bumh6 + _bumh7 
				* pow(tmgrid, 2)) * singrid + (0.60 * _bumh6 * (solargrid - 14930) / 3380.4));
}

// ***************************ecoregion 10020110 ***************************************** 
//vgtmp 2=lbalp, vgtmp3=lbpkl, vgtmp4=lbmhz, vgtmp5=lbsfz, vgtmp6=lbwhz, vgtmp7=lbsafz, 
//vgtmp8=lbgfz, vgtmp9=lbdf, vgtmp10=lbppz, vgtmp11=ubmh, vgtmp12=ubgf

// int v2ecoregion = con (((elevtmp < vgtmp10) and (elevtmp < vgtmp9)), 5, 33); // NonForest
int v3ecoregion = con (((elevtmp < vgtmp9) and (elevtmp > vgtmp10)), 10, 33); // PP
int v4ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp6) and (elevtmp < vgtmp7)), 14, 33); // DF
int v7ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp5) and (elevtmp < vgtmp7)), 19, 33); // WH
int v12ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp5) and (elevtmp < vgtmp11)), 19, 33); // WH
int v8ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp11) and (elevtmp < vgtmp3)), 22, 33); // PSF
int v9ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
int v10ecoregion = con (elevtmp < vgtmp3, 25, 33); // SAF
int v11ecoregion = con (elevtmp < vgtmp2, 32, 33); // PK

int vz_ecoregion = v3ecoregion;
if (v4ecoregion<vz_ecoregion) vz_ecoregion = v4ecoregion;
if (v7ecoregion<vz_ecoregion) vz_ecoregion = v7ecoregion;
if (v8ecoregion<vz_ecoregion) vz_ecoregion = v8ecoregion;
if (v9ecoregion<vz_ecoregion) vz_ecoregion = v9ecoregion;
if (v10ecoregion<vz_ecoregion) vz_ecoregion = v10ecoregion;
if (v11ecoregion<vz_ecoregion) vz_ecoregion = v11ecoregion;
if (v12ecoregion<vz_ecoregion) vz_ecoregion = v12ecoregion;

return((PNVcode)vz_ecoregion);
} // end of Subregion10020110()
*/

PNVcode PNVmodel::Subregion10020202()
{ // Wenatchee National Forest
	//*******************************************************************
	//*******************************************************************
	//   Alpine (ALPZ) Lower Boundary
	//*******************************************************************
	{
		float _blalp1 = 400;
		float _blalp2 = 145;
		float _blalp3 = -20;
		float _blalp4 = -16;
		float _blalp5 = 250;
		float _blalp6 = -2;
		//cad
		float _blalp11 = -0.3;
		//fog
		float _blalp12 = -500;

		vgtmp2 = int(((( _blalp11 * cadgrid + _blalp12 * foggrid + 
							_blalp1 + _blalp2 * msltgrid) + 
						(_blalp3 * tmgrid)) + (_blalp4 * prcpgrid)) + 
				(0.40 * _blalp5 + _blalp6 * prcpgrid) * singrid + 
				(0.60 * _blalp5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Parkland (PKZ) Lower Boundary
	//*******************************************************************
	{
		float _blpk1 = -800;
		float _blpk2 = 148;
		float _blpk3 = -50;
		float _blpk4 = -17;
		float _blpk5 = 180;
		float _blpk6 = 3;
		//cad
		float _blpk11 = -0.3;
		//fog
		float _blpk12 = -900;

		vgtmp3 = int((((  _blpk11 * cadgrid + _blpk12 * foggrid + 
							_blpk1 + _blpk2 * msltgrid) + (_blpk3 * tmgrid)) + 
					(_blpk4 * prcpgrid)) + 
				(0.40 * _blpk5 + _blpk6 * prcpgrid) * singrid + 
				(0.60 * _blpk5 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Mountain Hemlock (MH) Lower Boundary
	//*******************************************************************
	{
		float _blmh1 = -2300;
		float _blmh2 = 190;
		float _blmh3 = -180;
		float _blmh4 = -22;
		float _blmh5 = -2.5;
		float _blmh6 = 550;
		float _blmh7 = -30;
		float _blmh8 = -1.5;
		float _blmh9 = -0.1;
		//float _blmh11 -0.5
		//cad
		float _blmh11 = -1.0;
		//fog
		float _blmh12 = -450;

		vgtmp4 = int((((_blmh11 * cadgrid + _blmh12 * foggrid +  
							_blmh1 + _blmh2 * msltgrid) + (_blmh3 * pow(tmgrid, 1))) + 
					(_blmh4 + _blmh5 * pow(tmgrid, 1)) * prcpgrid) + 
				(( _blmh6 + _blmh7 * tmgrid) + 
				 (0.35 * _blmh8 + _blmh9 * tmgrid) * prcpgrid) * singrid + 
				(0.65 * _blmh8 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Pacific Silver Fir (PSF) Lower Boundary
	//********************************************************************
	{
		float _blpsf1 = -5650;
		float _blpsf2 = 245;
		float _blpsf3 = -250;
		float _blpsf4 = -27;
		float _blpsf5 = -2.6;
		float _blpsf6 = 600;
		float _blpsf7 = -20;
		float _blpsf8 = 15;
		float _blpsf9 = -2.2;
		float _blpsf11 = -1.2;
		float _blpsf12 = -900;

		vgtmp5 = int((((_blpsf11 * cadgrid + _blpsf12 * foggrid + 
							_blpsf1 + _blpsf2 * msltgrid) + (_blpsf3 * tmgrid)) + 
					(_blpsf4 + _blpsf5 * tmgrid) * prcpgrid) + 
				((0.30 * _blpsf6 + _blpsf7 * tmgrid) + 
				 (_blpsf8 + _blpsf9 * tmgrid) * prcpgrid) * singrid + 
				(0.70 * _blpsf8 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//    Western Hemlock (WH) Lower Boundary/Upper Grand Fir
	//*******************************************************************
	{
		float _blwh1 = -13500;
		float _blwh2 = 400;
		float _blwh3 = -100;
		float _blwh4 = 15;
		float _blwh5 = -20;
		float _blwh6 = 260;
		float _blwh7 = 55;
		float _blwh8 = 120;
		float _blwh9 = -18;
		//cad
		float _blwh11 = -2.0;
		//fog
		float _blwh12 = -900;

		vgtmp6 = int((((_blwh11 * cadgrid + _blwh12 * foggrid + 
							_blwh1 + _blwh2 * msltgrid) + (_blwh3 * tmgrid)) + 
					(_blwh4 + _blwh5 * tmgrid) * prcpgrid) + 
				((0.40 * _blwh6 + _blwh7 * tmgrid) + (_blwh8 + 
					_blwh9 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blwh6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//*******************************************************************
	//    Upper Grand Fir
	//*******************************************************************
	{
		float _bugf1 = -4100;
		float _bugf2 = 200;
		float _bugf3 = -140;
		float _bugf4 = -17;
		float _bugf5 = -5.2;
		float _bugf6 = 600;
		float _bugf7 = 90;
		float _bugf8 = 35;
		float _bugf9 = -10;
		//cad
		float _bugf11 = -1.0;
		//fog
		float _bugf12 = -1500;

		vgtmp12 = int((((_bugf11 * cadgrid + _bugf12 * foggrid + 
							_bugf1 + _bugf2 * msltgrid) + (_bugf3 * tmgrid)) + 
					(_bugf4 + _bugf5 * tmgrid) * prcpgrid) + 
				((0.40 * _bugf6 + _bugf7 * tmgrid) + (_bugf8 + _bugf9 * 
					tmgrid) * prcpgrid) * singrid + 
				(0.60 * _bugf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Subalpine Fir (SAFZ) Lower Boundary
	//*******************************************************************
	{
		float _blsaf1 = 2000;
		float _blsaf2 = 77;
		float _blsaf3 = -250;
		float _blsaf4 = 15;
		float _blsaf5 = -5;
		float _blsaf6 = 500;
		float _blsaf7 = -25;
		float _blsaf8 = -2.5;
		float _blsaf9 = 2;
		//cad
		float _blsaf11 = -0.4;
		//fog
		float _blsaf12 = 60;

		vgtmp7 = int((((   _blsaf11 * cadgrid + _blsaf12 * foggrid + 
							_blsaf1 + _blsaf2 * msltgrid) + (_blsaf3 * pow(tmgrid, .5))) + 
					(_blsaf4 + _blsaf5 * pow(tmgrid, .5)) * prcpgrid) + 
				(0.40 * _blsaf6 + _blsaf7 * tmgrid + (_blsaf8 + _blsaf9 * 
								      tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blsaf6 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//    Grand Fir (GF) Lower Boundary Grand Fir
	//*******************************************************************
	{
		float _blgf3 = 69.52;
		float _blgf2 = -7706.3;
		float _blgf1 = 220200;
		float _blgf4 = -365;
		float _blgf5 = 200;
		float _blgf6 = -132.85;
		float _blgf7 = 700;
		float _blgf8 = 25;
		float _blgf9 = 75;
		float _blgf10 = -17;
		float _blgf15 = 1.104;
		//cad
		float _blgf11 = -1.2;
		//fog
		float _blgf12 = -50;

		vgtmp8 = int(((((_blgf11 * cadgrid + _blgf12 * foggrid + 
								_blgf3 * pow(msltgrid, 2)) + _blgf2 * msltgrid + _blgf1) 
						+ _blgf4 * tmgrid) + (_blgf5 + _blgf6 * prcpgrid + _blgf15 * pow(prcpgrid,2))) + 
				((0.40 * _blgf7 + _blgf8 * tmgrid) + 
				 (_blgf9 + _blgf10 * tmgrid) * prcpgrid) * singrid + 
				(0.60 * _blgf7 * (solargrid - 14930) / 3380.4));
	}

	//*******************************************************************
	//*******************************************************************
	//   Douglas-fir (DF) Lower Boundary
	//*******************************************************************
	{
		float _bldf1 = 254500;
		float _bldf2 = -9110;
		float _bldf3 = 81.6;
		float _bldf4 = -120;
		float _bldf5 = 4500;
		float _bldf6 = -148.2;
		float _bldf7 = 200;
		float _bldf8 = 5;
		float _bldf9 = 300;
		float _bldf10 = 2500;
		float _bldf15 = 1.0724;
		//cad
		float _bldf11 = -0.9;
		//fog
		float _bldf12 = -500;

		vgtmp9 = int((_bldf11 * cadgrid + _bldf12 * foggrid + 
					_bldf1 + _bldf2 * msltgrid + _bldf3 * pow(msltgrid, 2) + 
					_bldf4 * tmgrid) + (_bldf5 + _bldf6 * prcpgrid + _bldf15 * pow(prcpgrid, 2)) + 
				(0.60 * _bldf7 + _bldf8 * pow(tmgrid, 0.5) + (_bldf9 + 
									      _bldf10 * tmgrid) / (prcpgrid + 10)) * singrid + 
				(0.40 * _bldf7 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Ponderosa Pine (PP) Lower Boundary
	//********************************************************************
	{
		float _blpp1 = 194800;
		float _blpp2 = -6917;
		float _blpp3 = 61.53;
		float _blpp4 = -250;
		float _blpp5 = 3888;
		float _blpp6 = -118.18;
		float _blpp7 = 280;
		float _blpp8 = 19;
		float _blpp9 = 30;
		float _blpp10 = -2.0;
		float _blpp15 = 0.743;
		//cad
		float _blpp11 = -0.8;
		//fog
		float _blpp12 = -300;

		vgtmp10 = int(( _blpp11 * cadgrid + _blpp12 * foggrid + 
					_blpp1 + _blpp2 * msltgrid + _blpp3 * pow(msltgrid, 2) + 
					_blpp4 * tmgrid + 
					_blpp5 + _blpp6 * prcpgrid + _blpp15 * pow(prcpgrid, 2)) + 
				((0.60 * _blpp7 + _blpp8 * tmgrid) + (_blpp9 + _blpp10 * tmgrid) * prcpgrid) * singrid + 
				(0.40 * _blpp7 * (solargrid - 14930) / 3380.4));
	}

	//********************************************************************
	//********************************************************************
	//    Mountain Hemlock (MH) Upper Boundary
	//********************************************************************
	{
		float _bumh1 = 16750;
		float _bumh2 = -190;
		float _bumh3 = -10;
		float _bumh4 = -77000;
		float _bumh5 = 155;
		float _bumh6 = -800;
		float _bumh7 = 15.0;
		//float _bumh8 = 2;
		//cad
		float _bumh11 = 1.3;
		//fog
		float _bumh12 = 600;

		vgtmp11 = int((((  _bumh11 * cadgrid + _bumh12 * foggrid + 
							_bumh1 + _bumh2 * msltgrid) + ( _bumh3 * tmgrid)) + ( _bumh4 + 
							_bumh5 * pow(tmgrid, 2)) / prcpgrid) + 
				(0.40 * _bumh6 + _bumh7 * pow(tmgrid, 2)) * singrid + 
				(0.60 * _bumh6 * (solargrid - 14930) / 3380.4));
	}

	//******************************************************************** 
	//vgtmp 2=lbalp, vgtmp3=lbpkl, vgtmp4=lbmhz, vgtmp5=lbsfz, vgtmp6=lbwhz, vgtmp7=lbsafz, 
	//vgtmp8=lbgfz, vgtmp9=lbdf, vgtmp10=lbppz, vgtmp11=ubmh, vgtmp12=ubgf

	// int v2ecoregion = con (((elevtmp < vgtmp10) and (elevtmp < vgtmp9)), 5, 33); // NonForest
	int v3ecoregion = con (((elevtmp < vgtmp9) and (elevtmp > vgtmp10) and (elevtmp < vgtmp6)), 10, 33); // PP
	int v4ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp6) and (elevtmp < vgtmp8) and (elevtmp < vgtmp7)), 14, 33); // DF
	int v5ecoregion = con (((elevtmp > vgtmp8) and (elevtmp < vgtmp6) and (elevtmp < vgtmp12) and (elevtmp < vgtmp3)), 16, 33); // GF
	int v7ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp6) and (elevtmp < vgtmp5) and (elevtmp < vgtmp7)), 19, 33); // WH
	int v8ecoregion = con (((elevtmp < vgtmp4) and (elevtmp < vgtmp7) and (elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 22, 33); // PSF
	int v9ecoregion = con (((elevtmp < vgtmp3) and (elevtmp < vgtmp11)), 23, 33); // MH
	int v10ecoregion = con (((elevtmp < vgtmp3) and (elevtmp > vgtmp11)), 25, 33); // SAF
	int v11ecoregion = con (((elevtmp < vgtmp2) and (elevtmp > vgtmp3) and (elevtmp > 3500)), 32, 33); // PK

	// int vz_ecoregion2 = v2ecoregion;
	// if (v3ecoregion<vz_ecoregion2) vz_ecoregion2 = v3ecoregion;
	int vz_ecoregion2 = v3ecoregion;
	if (v4ecoregion<vz_ecoregion2) vz_ecoregion2 = v4ecoregion;
	if (v5ecoregion<vz_ecoregion2) vz_ecoregion2 = v5ecoregion;
	if (v7ecoregion<vz_ecoregion2) vz_ecoregion2 = v7ecoregion;
	if (v8ecoregion<vz_ecoregion2) vz_ecoregion2 = v8ecoregion;
	if (v9ecoregion<vz_ecoregion2) vz_ecoregion2 = v9ecoregion;
	if (v10ecoregion<vz_ecoregion2) vz_ecoregion2 = v10ecoregion;
	if (v11ecoregion<vz_ecoregion2) vz_ecoregion2 = v11ecoregion;

	//---Renumber WH and DF to match numbers in table and shadesets---------
	int vz_ecoregion = con((vz_ecoregion2 == 17), 19, (con((vz_ecoregion2 == 18), 14, vz_ecoregion2)));

	return((PNVcode)vz_ecoregion);

} // end of Subregion10020202()
