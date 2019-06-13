package com.github.theLongLab.PRESM;

import java.util.Random;
import java.util.Vector;
import java.io.*;
import java.util.*;
import java.text.*;
import java.lang.*;

public class fasta {
	public Vector <String> ChrVec = new Vector<String>();
	public Vector <String> ChrSeqVec = new Vector<String>();
	public Vector <String> HomoeChrVec = new Vector<String>();
	public Vector <String> HomoChrSeqVec = new Vector<String>();

    public String [] description;
    public String [] sequence;


    public Vector <chrFasta> chrFastaVec = new Vector<chrFasta>();



	public fasta(String faFile) throws IOException{
		List desc= new ArrayList();
		List seq = new ArrayList();
		String temps="";

		try{
	        BufferedReader in     = new BufferedReader( new FileReader( faFile ) );
	        StringBuffer   buffer = new StringBuffer();
	        String         line   = in.readLine();

	        if( line == null )
	            throw new IOException( faFile + " is an empty file" );

	        if( line.charAt( 0 ) != '>' )
	            throw new IOException( "First line of " + faFile + " should start with '>'" );
	        else
	            desc.add(line);
	        for( line = in.readLine().trim(); line != null; line = in.readLine() )
		{
	            if( line.length()>0 && line.charAt( 0 ) == '>' )
		    {

			seq.add(buffer.toString());
			buffer = new StringBuffer();
			String chr=line;
			String[] chrv = chr.split(" ");
			desc.add(chrv[0]);
		    } else  {

	            	buffer.append( line.trim() );
		    	}
	        }
	        if( buffer.length() != 0 ){
	        	seq.add(buffer.toString());
	        }
	      }catch(IOException e)
	      {
	        System.out.println("Error when reading "+faFile);
	        e.printStackTrace();
	      }

		description = new String[desc.size()];
		sequence = new String[seq.size()];
		for (int i=0; i< seq.size(); i++)
		{
		  description[i]=(String) desc.get(i);
		  sequence[i]=(String) seq.get(i);
		}

		for (int i=0; i< seq.size()  ;i++){
			String[] splitedInfor = description[i].split(" ");
			String temp= splitedInfor[0];
			if (temp.length()!=0){
				if (temp.substring(0, 1).equals(">")) temp=temp.substring(1, temp.length());
			}
			ChrVec.add(temp );
			ChrSeqVec.add(sequence[i]);
		}

		for (int i=0;i<ChrVec.size();i++  ){
			chrFasta vc= new chrFasta();
			vc.assignment(ChrVec.get(i),ChrSeqVec.get(i));
			chrFastaVec.add(vc);
		}
	}

	public void ViewChr(String chr, int start, int end ) throws IOException{
		for (int i=0;i< ChrVec.size();i++){
			if (ChrVec.get(i).equals(chr) ){
				if (end > ChrSeqVec.get(i).length()) end = ChrSeqVec.get(i).length();
				System.out.println( "chromsome: "+chr+" From: "+start+" to: "+end );
				System.out.println( ChrSeqVec.get(i).substring(start-1, end) );
				return;
			}
		}
		System.out.println("Can not find chromsome: "+chr+" in the reference genome!");
	}

//MakeNewFastq

	public void MakeNewFastq( Vector <VcfChr> vec ) throws IOException{

		int index =0;
		String qual=         "888888888888888888888888888888888888888888888888888888888888888888888888888";
		for (int i=0; i< vec.size();i++){
			String chr = vec.get(i).ChrName;
			for (int j=0;j< ChrVec.size();j++){
				if (ChrVec.get(j).equals(chr) ){
					for (int k=0;k< vec.get(i).StartPositionVec.size();k++){
						for (int s=0;s< 70;s++){
							String seq= ChrSeqVec.get(j).substring( (vec.get(i).StartPositionVec.get(k)-70+s),
									(vec.get(i).StartPositionVec.get(k)+5+s)  );
							index++;
							String name= "@"+chr+"_"+ String.valueOf( vec.get(i).StartPositionVec.get(k) )+"_"+
									String.valueOf(vec.get(i).StartPositionVec.get(k)-69+s) +"_" + String.valueOf(index);
							System.out.println( name);
							System.out.println( seq);
							System.out.println( "+");
							System.out.println( qual);
						}
					}
				}
			}
		}
	}

	public void RandomGenerateMuation(int NumberOfMut, String type, int len, String OutputFile ) throws IOException{
		long  totalnumberofbase = 0;
		for (int i=0 ;i<chrFastaVec.size();i++ ){
			totalnumberofbase+= chrFastaVec.get(i).chrSeq.length();
//			System.out.println( totalnumberofbase );
		}
		long min = 0;
		long max = totalnumberofbase-1;
		Vector <Long> SedVec = new Vector<Long>();
		while (SedVec.size() < NumberOfMut) {
			Random random= new Random();
			int temp = (int)((max-200) /100);
			boolean flag = true;
			int hundred = random.nextInt(temp);
			int unit = random.nextInt(99);
			long longhundred =(long)hundred;
			long longunit =(long)unit;
			long num = 100*longhundred + longunit;
			for (int j=0; j< SedVec.size();j++){
				if (((SedVec.get(j) - num)<3)  &&  ((SedVec.get(j) - num)>-3))  {
					flag=false;
				}
			}
			if (flag ){
				SedVec.add(num);
			}
		}
		for (int i=0;i< SedVec.size();i++){
			for (int j=i;j< SedVec.size();j++){
				if (SedVec.get(i)> SedVec.get(j)){
					long temp = SedVec.get(i);
					SedVec.set(i, SedVec.get(j));
					SedVec.set(j, temp);
				}
			}
		}
		for (int i=0;i< SedVec.size();i++){
			System.out.println(  SedVec.get(i)  );
		}

		if (type=="snp"){
			for (int i=0;i<SedVec.size();i++ ){
				long num = SedVec.get(i);
				long start =0;
				long end =0;
				for (int j=0;j< chrFastaVec.size();j++){
					end += chrFastaVec.get(i).chrSeq.length();
					if (( num >=start) && (num<= end)){

					}
				}

			}
		}

//		while (SedVec.size() < NumberOfMut ) {
//			Random random= new Random();
//			boolean flag = false;
//			long s = random.nextInt(max)%(max-min+1)+min;
//			s= random.
//		}


	}


	private int abs(long l) {
		// TODO Auto-generated method stub
		return 0;
	}



}