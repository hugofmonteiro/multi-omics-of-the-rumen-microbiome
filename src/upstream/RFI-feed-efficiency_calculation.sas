*********************************************************************************	
*  Fábio Rumen Microbiome - Residual Feed Intake Calculations: May 2021			*
*********************************************************************************;
PROC IMPORT OUT= Work.RFI
            DATAFILE= "C:\Feed Efficiency and Rumen Microbiome\
Rumen Microbiome - USDA\Stats\SAS PROC & Files\1.RFIStats.csv"
            DBMS=CSV REPLACE;
     		GETNAMES=YES;
RUN;
DATA RFINew;
Set RFI;
IF ID = 'ID' THEN DELETE;
RUN;

PROC SORT DATA=RFINew;
by Seq;
PROC PRINT DATA=RFINew;
RUN;

*************************************************
*			Model for Residual DMI 				*
*************************************************;
* ID Seq FDAT Site Exp ExpN TRT TRTN Lact Parity FirstDIM LastDIM DOBS 
DMI DietNEL NELI Milk PCTF PCTP PCTL FatY ProtY LactY ECM NESec BW MBW BWC BCS BEC;

* Model for RFI: DMI = NESec MBW BEC Parity ;
*Random TRTN;

ODS PDF File='C:\Rumen Microbiome - USDA\
Stats\Outputs\RFI\1. RFI 40 to 120.PDF';

ODS HTML;
ODS Graphics On;
Title 'Fábio Rumen Microbiome - MIXED - Residual DM Intake from 40 to 120 DIM, kg/d';
RUN;
PROC MIXED  Data=RFINew COVTEST RATIO ITDETAILS MAXITER=500;
Class ID Seq EXP EXPN TRT TRTN Parity; 
MODEL DMI = NESec MBW BEC Parity / DDFM=KR S Residual Outp=RDMI;
Random TRTN /Solution;
Lsmeans Parity / PDIFF ADJUST=Tukey ADJDFE=Row;
RUN;

ODS Graphics Off;
ODS PDF Close;
QUIT;

PROC SORT DATA=RDMI; By Seq; 
PROC PRINT DATA=RDMI;
RUN;

