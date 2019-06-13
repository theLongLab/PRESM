package com.github.theLongLab.PRESM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import htsjdk.samtools.*;

public class Bam {


	public int CompareTwoStrings(String s1, String s2){
		int n=0;
		if (s1.length()!=s2.length()) return 100;
		for (int i=0;i< s1.length();i++){
			if (!s1.substring(i, i+1).equals(s2.substring(i, i+1))) n++;
		}
		return n;
	}

	public int AddNumberInString(String s){
		s=s+"*";
		int num=0;
		for (int i=0;i < s.length();i++){
			char si=s.charAt(i);
			if (Character.isDigit(si) ){
				for ( int j=i;j< s.length();j++){
					char sj=s.charAt(j);
					if (!Character.isDigit(sj) ){
						num+= Integer.parseInt( s.substring(i ,j) );
						i=j-1;
						break;
					}
				}
			}
		}
		return num;
	}

	//String[] splited = ContigLine.split("\t");

	public int ComputeInsSize(String s){
		int num=0;
		s="*"+s;
		for (int i=1; i<s.length();i++  ){
			if (s.substring(i, i+1).equals("I") ){
				for (int j=(i-1);j>-1;j--){
					char sj=s.charAt(j);
					if (!Character.isDigit(sj) ){
						num+= Integer.parseInt( s.substring(j+1,i) );
						break;
					}
				}
			}
		}
		return num;
	}


	public int[] ComputeMismatchAndDel(String s){

			int [] arr=new int[2];
			arr[0]=0;
			arr[1]=0;
			int num=0;
			for (int i=0;i<s.length();i++){
				char si=s.charAt(i);
				if ((si>='a' && si<='z') || (si>='A' && si<='Z')){
					num++;
				}
			}
			int DelNum=0;
			s=s+"9";
			for (int i=0;i< s.length();i++){
				if (s.subSequence(i, i+1).equals("^")){
					for ( int j=i+1;j< s.length();j++){
						char sj=s.charAt(j);
						if (!((sj>='a' && sj<='z') || (sj>='A' && sj<='Z'))){
							DelNum+= j-i-1 ;
							i=j-1;
							break;
						}
					}
				}
			}

			arr[0]=num-DelNum;
			arr[1]=DelNum;
			return arr;
	}

	public void SamSelect(String bFile, String outbFile) throws IOException{
		FileReader fr =new FileReader(bFile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		FileWriter mydata = new FileWriter(outbFile,true);
		PrintWriter fawrite = new PrintWriter(mydata);
		String tmp_name=null;
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("@")){
					int NumberOfSnp=10000;
					String[] splited = line.split("\t");
					String name= splited[0];
					String read=null;
					if (name.equals(tmp_name)){
						read= "_2";
					}else{
						tmp_name =name;
						read= "_1";
					}
					if (splited[2] .equals("*")){
						fawrite.write(name+read+"\t"+NumberOfSnp+"\n");
					}else{
						String str = splited[12];
						if (str.substring(0, 5).equals("MD:Z:")){
							String MDtag = str.substring(5, str.length());
			        		int MatchSize= AddNumberInString( MDtag) ;
			        		int MismatchSize = ComputeMismatchAndDel(MDtag)[0];
			        		int DelSize =  ComputeMismatchAndDel(MDtag)[1];
			        		int InsSize = ComputeInsSize( splited[5]  );
			        		NumberOfSnp= MismatchSize+ 4* DelSize+4* InsSize;
			        		fawrite.write(name+read+"\t"+NumberOfSnp+"\n");
						}
					}
				}
			}
		}
		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}


	public void FalseMap(String bFile, String outbFile) throws IOException{
		FileReader fr =new FileReader(bFile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		FileWriter mydata = new FileWriter(outbFile,true);
		PrintWriter fawrite = new PrintWriter(mydata);
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("@")){
					int NumberOfSnp=10000;
					String[] splited = line.split("\t");
					if (line.substring(0, 4) .equals("rand")){
//						fawrite.write(NumberOfSnp+"\n");
						NumberOfSnp=99;
					}else{
						String str= splited[0];
						String[] splited_str = str.split("_");
						String simu_chr = splited_str[0];
						int simu_pos= Integer.parseInt(  splited_str[1].toString());
						String map_chr= splited[2];
						int map_pos = Integer.parseInt(  splited[3].toString());
						if (( !simu_chr.equals(map_chr)) ||   ( (Math.abs (  simu_pos- map_pos ))  > 10000000 ) ){
							fawrite.write( line+"\n") ;
						}
					}
//		    		System.out.println( line) ;1
//		    		System.out.println(NumberOfSnp);
				}
			}
		}
		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}

	public void RevertSam(String bFile, String outbFile, String listfile) throws IOException{

		Vector <String> listVec = new Vector<String>();
		FileReader fr =new FileReader(listfile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					listVec.add(line);
				}
			}
		}
		bufferedreader.close();

		File bamFile = new File( bFile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        String file =outbFile;
		File writebamFile = new File (file);
		SAMFileWriter outputBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),
		true, writebamFile);
        SAMRecordIterator r = sr.iterator();

        while(r.hasNext()) {
        	SAMRecord rl= r.next();

        	boolean ok =true;
        	for (int i=0;i< listVec.size() ;i++ ){

        		if (rl.getReadName().equals( listVec.get(i))) {
        			ok=false;
        			System.out.println( rl.getReadName()   );
        		}
         	}
        	if (ok){
        		outputBam.addAlignment(rl );
        	}
        }
        r.close();
        sr.close();
        outputBam.close();

    }







	public void BamSelect(String bFile, String outbFile, int n) throws IOException{
		File bamFile = new File( bFile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        String file =outbFile;
		File writebamFile = new File (file);
		SAMFileWriter outputBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),
		true, writebamFile);
        SAMRecordIterator r = sr.iterator();
        while(r.hasNext()) {
        	SAMRecord rl= r.next();
        	int MatchSize=0;
        	int MismatchSize=0;
        	int InsSize=0;
        	int DelSize=0;
        	int NumberOfSnp=0;// del=4 snp, invertion =4 snp\

        	if (rl.getAttribute("MD")==null  ) {
        		NumberOfSnp+=1000;
        	} else {
        		String MDtag=rl.getAttribute("MD").toString();
        		MatchSize= AddNumberInString( MDtag) ;
        		MismatchSize = ComputeMismatchAndDel(MDtag)[0];
        		DelSize =  ComputeMismatchAndDel(MDtag)[1];
        		InsSize = ComputeInsSize(rl.getCigarString() );
        	}
        	NumberOfSnp+= MismatchSize+ 4* DelSize+4* InsSize;

        	if  ((  rl.getSAMString().toString().length()!=0)  &&   (rl.getSAMString().toString().substring(0, 1).equals("@"))){
        		NumberOfSnp=0;
//        		if (rl.getAttribute("MD")==null  )  {
//        			System.out.println( rl.getSAMString()) ;
//        		}

        	}
//    		if ((rl.getAttribute("MD")==null  ) &&  (!(rl.getSAMString().toString().substring(0, 1).equals("@"))))   {
//    			if (!rl.getCigarString().equals("*") )
//    				System.out.println( rl.getSAMString()) ;
//    		}
//        	System.out.println( rl.getSAMString()) ;
//        	System.out.println( rl.getCigar().numCigarElements()) ;
//        	System.out.println( rl.getCigar().getCigarElements() ) ;

/*    		FileWriter mydata = new FileWriter(fileTMP,true);
    		PrintWriter fawrite = new PrintWriter(mydata);*/

        	if (rl.getCigarString().equals("*") ) NumberOfSnp=100;
        	boolean ok=false;
        	if (NumberOfSnp<(n+1)  ){
//        	if (NumberOfSnp!=-9999  ){
        		outputBam.addAlignment(rl );
        		ok=true;
        	}else{

//        		System.out.println( rl.getSAMString()) ;
        	}


        }
        r.close();
        sr.close();
        outputBam.close();
//        final File inputSamOrBamFile;
//        final File outputSamOrBamFile;
//        final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
//        true, outputSamOrBamFile);
//        for (final SAMRecord samRecord : inputSamOrBamFile) {// Convert read name to upper case.
//        	samRecord.setReadName(samRecord.getReadName().toUpperCase());
//        	outputSam.addAlignment(samRecord);
//        	}
//        	outputSam.close();
//        	inputSam.close();
//		}
    }

	public void SplitBam(String bFile, String outbFile, int n, String library) throws IOException{
		File bamFile = new File( bFile);
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        String file =outbFile;
		File writebamFile = new File (file);
		SAMFileWriter outputBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),
		true, writebamFile);
        SAMRecordIterator r = sr.iterator();
        Vector <String> LibStrVec = new Vector<String>();
        SAMRecord rl= r.next();


        String header = rl.getHeader().getSAMString() ;
        String[] HeaderSplit = header.split("\n");
        for (int i=0;i<HeaderSplit.length;i++ ){
        	String ss = HeaderSplit[i];
        	if (ss.length()>3){
        		if (ss.substring(0, 3)=="@RG"){
        			String[] lineSplit = ss.split("\t");
        			String id ="";
        			String lib ="";

        		}
        	}
		}
//        while(r.hasNext()) {
//        	SAMRecord rl= r.next();
////        	System.out.println(rl.getHeader().getSAMString() );
//
//        }
	}


}

/*
public void convertReadNamesToUpperCase(
final File inputSamOrBamFile,
final File outputSamOrBamFile
)
{
final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
final SAMFileWriter outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
true, outputSamOrBamFile);

for (final SAMRecord samRecord : inputSam) {
// Convert read name to upper case.
samRecord.setReadName(samRecord.getReadName().toUpperCase());
outputSam.addAlignment(samRecord);
}
outputSam.close();
inputSam.close();
}
(...)
*/