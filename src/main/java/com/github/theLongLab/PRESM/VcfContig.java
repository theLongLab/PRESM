package com.github.theLongLab.PRESM;

public class VcfContig extends vcf{

	public VcfContig(String ContigLine, int flag){
//		System.out.println( ContigLine  );
		String[] splited = ContigLine.split("\t");

//		for (int i=0;i<  splited.length;i++){
//			System.out.println( splited[i]+" " );
//		}
		ChrName=splited[0];
		StartPoint=Integer.parseInt(splited[1]);
		SNPid=splited[2];
		RefStr=splited[3];
		EndPoint=StartPoint+RefStr.length()-1;

		AltStr=splited[4];
		int p =AltStr.indexOf(",");
		if (p!=-1) AltStr=AltStr.substring(0, p);
		Qual=splited[5];
		Filter= splited[6];
		Info = splited[7];

		SourceFlag=flag;

		IsPindelVcf=false;//*******For Vcf File Made by Pindel
		int len = splited.length;

		String gtfield= splited[len-1];
		String[] gtvec=  gtfield.split(":");
		String gt = gtvec[0];
		if ((gt.equals("0/1")) || (gt.equals("1/0")) || (gt.equals("1|0"))||(gt.equals("0|1")) ){
			Ishetero=true;
		}else{
			Ishetero=false;
		}
		gtfield= splited[len-2];
		String[] gtvec2=  gtfield.split(":");
		gt = gtvec2[0];
		if ((gt.equals("0/1")) || (gt.equals("1/0")) || (gt.equals("1|0"))||(gt.equals("0|1")) ){
			Ishetero=true;
		}
//		if (IsPindelVcf){
		if (flag==1){
//			if (!ChrName.substring(0, 1).equals("c"))	ChrName="chr"+ChrName;
//			String[] splitedInfor = Info.split(";");

/*
			for (int i=0;i<  splitedInfor.length;i++){
				SVLen= Integer.parseInt(splitedInfor[2].substring(6, splitedInfor[2].length()));
				NTLen=Integer.parseInt(splitedInfor[4].substring(6, splitedInfor[4].length()));
				SVType= splitedInfor[3].substring(7, splitedInfor[3].length()) ;
			}
*/

//SV Type:Copied from pindel2vcftest.cpp which located in src folder of Pindel
//outFile << "##ALT=<ID=DEL,Description=\"Deletion\">" << endl; /*EWL040311: probably not needed, as our calls are precise so we should rather give the exact replacing sequence instead of a label. */
//outFile << "##ALT=<ID=DUP,Description=\"Duplication\">" << endl;
//outFile << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">" << endl;
//outFile << "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl;
//outFile << "##ALT=<ID=INV,Description=\"Inversion\">" << endl;
//outFile << "##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl;


		}
	}

	public String ChrName;
	public int StartPoint;
	public int EndPoint;
	public int AltEndPoint;
	public String SNPid;
	public String Qual;
	public String Filter;
	public String RefStr;
	public String AltStr;
	public String Info;
	public boolean IsGatkVcf;
	public boolean IsPindelVcf;
	public int SVLen;
	public String SVType;
	public int NTLen;
	public int SourceFlag;
	public boolean Ishetero;
}
