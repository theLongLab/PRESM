package com.github.theLongLab.PRESM;


import java.lang.*;
import java.io.*;
import java.util.*;
import java.text.*;
import htsjdk.samtools.*;


public class vcf{

	public String VcfFilePath;
	public String GatkVcfFilePath;
	public String PindelVcfFilePath;
	public Vector <VcfContig> VcfConfigVec = new Vector<VcfContig>();
	public Hashtable<String,Integer> chr_index_hm  =new Hashtable<String,Integer>();
	public Vector <VcfChr> VcfChrVec = new Vector<VcfChr>();
	String TransferRead(String seq, String qual, String cigar, int pos, int[] startvec, int[] endvec, String[] refvec, String[] altvec){
		if (cigar.equals("*"))  return seq+"\t"+qual;
    	String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
    	int[] cigarnums = new int[cigarnumsTmp.length];
        for(int i = 0; i < cigarnumsTmp.length; i++){
        	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
        }
        String cigarlettersTmp[] = cigar.split("[0-9]+");//.toString().toCharArray();
        String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
        String[] basevec = new String[seq.length()];
        String[] qualvec= new String[seq.length()];
        Boolean[] insertarea= new Boolean[seq.length()];
        int[] indexvec = new int[seq.length()];
        for (int i=0;i<basevec.length;i++ ){
        	basevec[i]= seq.substring(i, i+1);
        	qualvec[i]= qual.substring(i,i+1);
        	indexvec[i]= pos+i;
        	insertarea[i]=false;
        }
        int position =pos;
        boolean hasunknownoperation=false;
        for (int i=0;i<cigarletters.length;i++ ){
        	if ((cigarletters[i].equals("S")) ||   (cigarletters[i].equals("S")) ) {
        		for (int j=position;j<(position+cigarnums[i]);j++ ){
        			if ((j -pos)< indexvec.length  ) {
        				indexvec[j-pos]=-1;
        			}
        		}
        		for (int j=(position+cigarnums[i]);j<(pos+indexvec.length);j++ ){
        				indexvec[j-pos]=indexvec[j-pos]-cigarnums[i] ;
        		}
        		position = position + cigarnums[i];
        	}
        	if (cigarletters[i].equals("M"))  {
        		position = position + cigarnums[i];
        	}
        	if (cigarletters[i].equals("D"))  {
        		for (int j=position; j< (pos+indexvec.length );j++){
        			indexvec[j-pos]=indexvec[j-pos]+cigarnums[i] ;
        		}
        	}
        	if (cigarletters[i].equals("I")) {
        		for (int j=position; j< (position+cigarnums[i]);j++){
        				insertarea[j-pos]=true;
        		}
        		position= position + cigarnums[i];
        		for (int j=position; j< (pos+indexvec.length);j++){
        			indexvec[j-pos]=indexvec[j-pos]-cigarnums[i] ;
        		}
        	}
        	if ((cigarletters[i].equals("N"))  || (cigarletters[i].equals("P"))   || (cigarletters[i].equals("H"))    ){
        		hasunknownoperation=true;
        	}
        }

        if (hasunknownoperation){
        	return seq+"\t"+qual;
        }

        boolean hasstar=false;
        for (int i=0;i<startvec.length;i++ ){
        	boolean check=true;
        	int firstindex=-1;
        	for (int j=0;j< refvec[i].length();j++){
        		int refbasepoi= startvec[i]+j;
        		String refbase= refvec[i].substring(j, j+1);
        		int num=0;
        		for (int s=0;s< indexvec.length;s++){
        			if (indexvec[s]==refbasepoi){
        				if (firstindex==-1) firstindex =s;
        				num++;
        				if (!basevec[s].equals(refbase) ){
        					check=false;
        				}
        				if  ((s+1) < insertarea.length){
        					if (insertarea[s+1]){
        						check=false;
        					}
        				}
        			}
        		}
        		if (num>1)  check=false;
        	}
        	if (firstindex==-1) check=false;
        	if (check){
//        		System.out.println( pos+" "+ seq+ " "+cigar);
//        		System.out.println(startvec[i]+" "+  refvec[i]+" "+altvec[i]);
        		for (int j=0;j<refvec[i].length();j++ ){
        			int temp_pos= startvec[i]+j;
        			for (int s=0;s< indexvec.length;s++){
        				if (indexvec[s]== temp_pos){
//        					System.out.println("BBC"+"\t"+temp_pos+"\t"+refvec[i]+"\t"+ refvec[i].substring(j, j+1)+"\t"+ altvec[i]+"\t"+ basevec[s]);
        				}
        			}
        		}

        		hasstar=true;
        		for (int j=0;j< refvec[i].length();j++){
            		int refbasepoi= startvec[i]+j;
//            		String refbase= refvec[i].substring(j, j+1);
            		for (int s=0;s< indexvec.length;s++){
            			if (indexvec[s]==refbasepoi){
            				if (s==firstindex)  {
            					basevec[s]=altvec[i];
            					String temp="";
            					for (int r=0;r< altvec[i].length();r++)
            						temp=temp+"I";
            					qualvec[s]=temp;
            				}else{
            					basevec[s]="*";
            					qualvec[s]="*";
            				}
            			}
            		}
        		}
        	}else{
        		for (int j=0;j<refvec[i].length();j++ ){
        			int temp_pos= startvec[i]+j;
        			for (int s=0;s< indexvec.length;s++){
        				if (indexvec[s]== temp_pos){
        				}
        			}
        		}
        	}
        }

        String convert_seq="";
        String convert_qual="";

        for (int i=0;i< basevec.length;i++){
        	if (!basevec[i].equals( "*")){
        		convert_seq= convert_seq+ basevec[i];
        		convert_qual= convert_qual+ qualvec[i];
        	}
        }
//        if (!hasstar) System.out.println("Not Match the Reference Genome #######################" );
        if ((convert_seq.length()< (2*seq.length()) )  && ( convert_seq.length()>18 ) )  {
//        	if (hasstar)System.out.println("Match the Reference Genome **********************" );
        	return convert_seq+"\t"+convert_qual;
        }else{
        	return  seq+"\t"+qual;
        }
	}

	void AlterVcf(String file, fasta ff) throws IOException{ //If just one vcf is provided
		String TumorNewVcfFile="/home/chencao/Desktop/vcf/new.vcf";
		FileWriter mydata = new FileWriter(TumorNewVcfFile,true);
		PrintWriter VcfWriter = new PrintWriter(mydata);


		FileReader fr2 =new FileReader(file);
		BufferedReader bufferedreader2= new BufferedReader(fr2);
		String line="";
		while ( (line =bufferedreader2.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#")){
					VcfWriter.write(line+"\n");
				}
				if (!line.substring(0, 1).equals("#")){

					boolean boring=false;
					String[] splited = line.split("\t");
					if (splited[4].equals("<DUP:TANDEM>")){
						String chr=splited[0];
						int start=Integer.parseInt(splited[1]);
						String[] infosplited = splited[7].split("=");
						int end= start + Integer.parseInt(infosplited[infosplited.length-1])-1;
						boring=true;
						for (int i=0;i< ff.ChrVec.size();i++){
							if (ff.ChrVec.get(i).equals(chr)){
								splited[4]=ff.ChrSeqVec.get(i).substring(start-1, end);
							}
						}
						for (int i=0;i< splited.length;i++){
							if (i!= (splited.length-1) ){
								VcfWriter.write(splited[i]+"\t");
							}else {
								VcfWriter.write(splited[i]+"\n");
							}
						}
					}
					if (splited[4].equals("<INV>")){
						String chr=splited[0];
						int start=Integer.parseInt(splited[1]);
						String[] infosplited = splited[7] .split("=");
						int end= start + Integer.parseInt(infosplited[infosplited.length-1])-1;
						boring=true;
						for (int i=0;i< ff.ChrVec.size();i++){
							if (ff.ChrVec.get(i).equals(chr)){
								splited[4]=ff.ChrSeqVec.get(i).substring(start-1, end);
								StringBuffer sb=new StringBuffer( splited[4]);
								sb.reverse();
								splited[3]=sb.toString();
							}
						}
						for (int i=0;i< splited.length;i++){
							if (i!= (splited.length-1) ){
								VcfWriter.write(splited[i]+"\t");
							}else {
								VcfWriter.write(splited[i]+"\n");
							}
						}
					}
					if (splited[4].indexOf(",")!=-1){
						String[] no4vec= splited[4].split(",");
						String[] no7vec= splited[7].split(",");
						splited[4]= no4vec[0];
						splited[7]= no7vec[0];
						splited[9]="1|1";
						for (int i=0;i< splited.length;i++){
							if (i!= (splited.length-1) ){
								VcfWriter.write(splited[i]+"\t");
							}else {
								VcfWriter.write(splited[i]+"\n");
							}
						}
						boring=true;
					}
					if (!boring)  VcfWriter.write(line+"\n");


				}
			}
		}
		VcfWriter.close();
	}

	void readVcfFromFile(String file) throws IOException{ //If just one vcf is provided
		VcfFilePath=file;
		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					VcfConfigVec.add(vc);
				}
			}
		}
		bufferedreader.close();
		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {
				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}
		for (int i=0; i< ChrVec.size();i++){
//			System.out.println(ChrVec.get(i));
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}
	}

	boolean meetcondition (VcfContig vc, String gt, Vector <String> ChrVec,Vector <Integer> StartVec,
			Vector <Integer> EndVec){
		boolean flag=false;
		for (int i=0;i< ChrVec.size();i++ ){
			if (ChrVec.get(i).equals(vc.ChrName)) {
				if ((vc.StartPoint>=StartVec.get(i)) &&  (vc.EndPoint<=EndVec.get(i))){
					flag=true;
				}
			}
		}
		if (!flag)return false;
		else if (gt ==null) return true;
		else if ((gt.equals("homo")) && (!vc.Ishetero)) return true;
		else if ((gt.equals("heter")) && (vc.Ishetero)) return true;
		return false;
	}

	void readVcfFromFile(String file, String gt, String interval, int len) throws IOException{
		String line=null;
		Vector <String> IntervalChrVec = new Vector<String>();
		Vector <Integer> IntervalStartVec = new Vector<Integer>();
		Vector <Integer> IntervalEndVec = new Vector<Integer>();
		if ( interval!=null){
			FileReader ir =new FileReader(interval);
			BufferedReader intervalreader= new BufferedReader(ir);
			while ( (line =intervalreader.readLine())!=null ){
				int min=1;
				int max= Integer.MAX_VALUE;
				String chr=null;
				if (0!=line.length()){
					String[] splited = line.split(":");
					if (splited.length==1){
						chr= splited[0];
					}else{
						String range = splited[1];
						chr= splited[0];
						String[] splited_range = range.split("-");
						if (splited_range.length==1){
							min = Integer.parseInt( splited_range[0]  );
						}else{
							min= Integer.parseInt( splited_range[0]  );
							max= Integer.parseInt( splited_range[1]  );
						}
					}
				}
				IntervalChrVec.add(chr);
				IntervalStartVec.add(min);
				IntervalEndVec.add(max);
			}
			intervalreader.close();
		}
		VcfFilePath=file;
		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					if (interval!=null){
						if (meetcondition(vc, gt, IntervalChrVec,IntervalStartVec, IntervalEndVec )){
							if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
								VcfConfigVec.add(vc);
						}
					}else{
						if (gt==null){
							if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
								VcfConfigVec.add(vc);
						}else if (gt.equals("homo")){
							if (!vc.Ishetero){
								if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
									VcfConfigVec.add(vc);
							}
						}else if (gt.equals("heter")){
							if (vc.Ishetero){
								if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
									VcfConfigVec.add(vc);
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {

				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}
		for (int i=0; i< ChrVec.size();i++){
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}

	}


	void readVcfFromFile_Meth(String file, String gt, String interval, int len) throws IOException{
		String line=null;
		Vector <String> IntervalChrVec = new Vector<String>();
		Vector <Integer> IntervalStartVec = new Vector<Integer>();
		Vector <Integer> IntervalEndVec = new Vector<Integer>();
		if ( interval!=null){
			FileReader ir =new FileReader(interval);
			BufferedReader intervalreader= new BufferedReader(ir);
			while ( (line =intervalreader.readLine())!=null ){
				int min=1;
				int max= Integer.MAX_VALUE;
				String chr=null;
				if (0!=line.length()){
					String[] splited = line.split(":");
					if (splited.length==1){
						chr= splited[0];
					}else{
						String range = splited[1];
						chr= splited[0];
						String[] splited_range = range.split("-");
						if (splited_range.length==1){
							min = Integer.parseInt( splited_range[0]  );
						}else{
							min= Integer.parseInt( splited_range[0]  );
							max= Integer.parseInt( splited_range[1]  );
						}
					}
				}
				IntervalChrVec.add(chr);
				IntervalStartVec.add(min);
				IntervalEndVec.add(max);
			}
			intervalreader.close();
		}
		VcfFilePath=file;
		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					if (interval!=null){
						if (meetcondition(vc, gt, IntervalChrVec,IntervalStartVec, IntervalEndVec )){
							if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
								VcfConfigVec.add(vc);
						}
					}else{
						if (gt==null){
							if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
								VcfConfigVec.add(vc);
						}else if (gt.equals("homo")){
							if (!vc.Ishetero){
								if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
									VcfConfigVec.add(vc);
							}
						}else if (gt.equals("heter")){
							if (vc.Ishetero){
								if ((vc.AltStr.length()<len ) && (vc.RefStr.length()<len ))
									VcfConfigVec.add(vc);
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {
				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}
		for (int i=0; i< ChrVec.size();i++){
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}
		for (int i=0; i< VcfChrVec.size();i++){
			chr_index_hm.put(VcfChrVec.get(i).ChrName, i);
		}
	}

	void readVcfFromFile_Meth(String file, String gt, String interval) throws IOException{
		String line=null;
		Vector <String> IntervalChrVec = new Vector<String>();
		Vector <Integer> IntervalStartVec = new Vector<Integer>();
		Vector <Integer> IntervalEndVec = new Vector<Integer>();
		if ( interval!=null){
			FileReader ir =new FileReader(interval);
			BufferedReader intervalreader= new BufferedReader(ir);
			while ( (line =intervalreader.readLine())!=null ){
				int min=1;
				int max= Integer.MAX_VALUE;
				String chr=null;
				if (0!=line.length()){
					String[] splited = line.split(":");
					if (splited.length==1){
						chr= splited[0];
					}else{
						String range = splited[1];
						chr= splited[0];
						String[] splited_range = range.split("-");
						if (splited_range.length==1){
							min = Integer.parseInt( splited_range[0]  );
						}else{
							min= Integer.parseInt( splited_range[0]  );
							max= Integer.parseInt( splited_range[1]  );
						}
					}
				}
				IntervalChrVec.add(chr);
				IntervalStartVec.add(min);
				IntervalEndVec.add(max);
			}
			intervalreader.close();
//			for (int i=0;i< IntervalChrVec.size() ;i++ ){
//				System.out.println(IntervalChrVec.get(i)+"\t"+ +IntervalStartVec.get(i)+"\t"+IntervalEndVec.get(i));
//			}
		}
		VcfFilePath=file;

		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					if (interval!=null){
						if (meetcondition(vc, gt, IntervalChrVec,IntervalStartVec, IntervalEndVec )){
							VcfConfigVec.add(vc);
						}
					}else{
						if (gt==null){
							VcfConfigVec.add(vc);
						}else if (gt.equals("homo")){
							if (!vc.Ishetero){
								VcfConfigVec.add(vc);
							}
						}else if (gt.equals("heter")){
							if (vc.Ishetero){
								VcfConfigVec.add(vc);
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
//		for (int i=0;i< VcfConfigVec.size() ;i++ ){
//			System.out.println(VcfConfigVec.get(i).ChrName+"\t"+ +VcfConfigVec.get(i).StartPoint+"\t"+VcfConfigVec.get(i).RefStr);
//		}
		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {
				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}
		for (int i=0; i< ChrVec.size();i++){
//			System.out.println(ChrVec.get(i));
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}
	}


	void readVcfFromFile(String file, String gt, String interval) throws IOException{
		String line=null;
		Vector <String> IntervalChrVec = new Vector<String>();
		Vector <Integer> IntervalStartVec = new Vector<Integer>();
		Vector <Integer> IntervalEndVec = new Vector<Integer>();
		if ( interval!=null){
			FileReader ir =new FileReader(interval);
			BufferedReader intervalreader= new BufferedReader(ir);
			while ( (line =intervalreader.readLine())!=null ){
				int min=1;
				int max= Integer.MAX_VALUE;
				String chr=null;
				if (0!=line.length()){
					String[] splited = line.split(":");
					if (splited.length==1){
						chr= splited[0];
					}else{
						String range = splited[1];
						chr= splited[0];
						String[] splited_range = range.split("-");
						if (splited_range.length==1){
							min = Integer.parseInt( splited_range[0]  );
						}else{
							min= Integer.parseInt( splited_range[0]  );
							max= Integer.parseInt( splited_range[1]  );
						}
					}
				}
				IntervalChrVec.add(chr);
				IntervalStartVec.add(min);
				IntervalEndVec.add(max);
			}
			intervalreader.close();
//			for (int i=0;i< IntervalChrVec.size() ;i++ ){
//				System.out.println(IntervalChrVec.get(i)+"\t"+ +IntervalStartVec.get(i)+"\t"+IntervalEndVec.get(i));
//			}
		}
		VcfFilePath=file;
		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					if (interval!=null){
						if (meetcondition(vc, gt, IntervalChrVec,IntervalStartVec, IntervalEndVec )){
							VcfConfigVec.add(vc);
						}
					}else{
						if (gt==null){
							VcfConfigVec.add(vc);
						}else if (gt.equals("homo")){
							if (!vc.Ishetero){
								VcfConfigVec.add(vc);
							}
						}else if (gt.equals("heter")){
							if (vc.Ishetero){
								VcfConfigVec.add(vc);
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
//		for (int i=0;i< VcfConfigVec.size() ;i++ ){
//			System.out.println(VcfConfigVec.get(i).ChrName+"\t"+ +VcfConfigVec.get(i).StartPoint+"\t"+VcfConfigVec.get(i).RefStr);
//		}
		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {
				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}
		for (int i=0; i< ChrVec.size();i++){
//			System.out.println(ChrVec.get(i));
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}
	}

	void readVcfFromFile(String file1, String file2) throws IOException{ //gatk vcf file for file1 and pindel vcf file for file2
		GatkVcfFilePath=file1;
		PindelVcfFilePath=file2;
		FileReader fr2 =new FileReader(file2);
		BufferedReader bufferedreader2= new BufferedReader(fr2);
		String line="";
		while ( (line =bufferedreader2.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,1); //1 for gatk
					VcfConfigVec.add(vc);
				}
			}
		}
// If the indel or snp located in the region of SV assigned by pindel, the indel or snp will be discarded
		FileReader fr =new FileReader(file1);
		BufferedReader bufferedreader= new BufferedReader(fr);
		line="";
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (!line.substring(0, 1).equals("#")){
					VcfContig vc= new VcfContig(line,0);
					VcfConfigVec.add(vc);
				}
			}
		}
		bufferedreader.close();

		bufferedreader2.close();

		Vector <String> ChrVec = new Vector<String>();
		for (int i=0;i< VcfConfigVec.size() ;i++ ){
			if  (!ChrVec.contains(VcfConfigVec.get(i).ChrName) ) {
				ChrVec.add( VcfConfigVec.get(i).ChrName);
			}
		}



		for (int i=0; i< ChrVec.size();i++){
			System.out.println(ChrVec.get(i));
			VcfChr vcr= new VcfChr(ChrVec.get(i),VcfConfigVec  );
			VcfChrVec.add(vcr);
		}
	}

	public int [] randomCommon(int min, int max, int n){
		if ( (n> (max-min+1) ) || (max<min) ){
			return null;
		}
		int [] result =new int [n];
		int count=0 ;
		while (count <n){
			int num = (int ) (Math.random() *(max-min+min));
			boolean flag=true;
			for (int j=0; j< n;j++){
				if (num == result[j]){
					flag=false;
					break;
				}
			}
			if (flag){
				result[count]=num;
				count++;
			}
		}
		return result;
	}


	void SelectVcf ( String dbsbpFile) throws IOException{
		FileReader fr =new FileReader(dbsbpFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		int count =0;
		int [] rand = randomCommon (56, 149120263, 50000);

		for (int i =0;i< rand.length;i++){
			for (int j=i;j< rand.length;j++){
				if (rand[i]>rand[j]){
					int temp= rand[i];
					rand[i]=rand[j];
					rand[j] =temp;
				}
			}
		}
		for (int i =0;i< rand.length;i++){
			System.out.println(rand[i]+"\t" );
		}

		String file="/home/chencao/Desktop/random.vcf";
		FileWriter vcfdata = new FileWriter(file,true);
		PrintWriter vcfwrite = new PrintWriter(vcfdata);


		while ( (line =bufferedreader.readLine())!=null ){
			count++;
			if (line.substring(0, 1).equals("#")){
				vcfwrite.write(line+"\n");
			}

			for (int i =0; i< rand.length;i++){
				if ( rand[i] == count)  {
					String[] splited = line.split("\t");
					if ((splited[3].length()==1 ) && (splited[4].length()==1 )){
						vcfwrite.write(line+"\n");
						break;
					}
				}
			}
		}

		bufferedreader.close();
		vcfwrite.flush();
		vcfwrite.close();
	}

	void MakePersonalizedVairants(Vector <VcfChr> vec,String TumorVcfFile,  String TumorNewVcfFile
			,boolean removeduplicates) throws IOException{
		FileWriter mydata = new FileWriter(TumorNewVcfFile,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		FileReader fr =new FileReader(TumorVcfFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		int LastPos=-1;
		int startpoint=0;
		String flagchr=null;
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#"))  {
					VcfWriter.write(line+"\n")  ;
				}else{
					String[] splitedInfor = line.split("\t");
					String chr = splitedInfor[0].substring(0, splitedInfor[0].length());
					int position=Integer.parseInt(splitedInfor[1]);
					int newpos=position;
					int move=0;
					if (!chr.equals(flagchr)){
						LastPos=-1;
						startpoint=0;
						flagchr=chr;
					}
					for (int i=0; i< vec.size();i++){
						if (chr.equals(vec.get(i).ChrName)){
							for (int j= startpoint; j< vec.get(i).AlterStartVec.size();j++){
									if (position > vec.get(i).AlterEndVec.get(j) ){
										move=vec.get(i).AlterDisVec.get(j)  ;
										startpoint=j;
									}else{
										break;
									}
							}
						}
					}
					newpos+=move;
					if  (((removeduplicates==true) && (newpos > LastPos ) )  || (!removeduplicates)){
						LastPos = newpos;
						VcfWriter.write(chr+"\t");
						VcfWriter.write(newpos+"\t");
						for (int i=2;i<splitedInfor.length;i++ ){
							if (i!= (splitedInfor.length-1) ){
								VcfWriter.write(splitedInfor[i]+"\t");
							}else {
								VcfWriter.write(splitedInfor[i]+"\n");
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		VcfWriter.close();
	}

	void MakeNewVcf(Vector <VcfChr> vec,String TumorVcfFile,  String TumorNewVcfFile) throws IOException{
		FileWriter mydata = new FileWriter(TumorNewVcfFile,true);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		FileReader fr =new FileReader(TumorVcfFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		int LastPos=-1;
		int startpoint=0;
		String flagchr=null;
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#"))  {
					VcfWriter.write(line+"\n")  ;
				}else{
					String[] splitedInfor = line.split("\t");
					String chr = splitedInfor[0].substring(0, splitedInfor[0].length());
					int position=Integer.parseInt(splitedInfor[1]);
					int newpos=position;
					int move=0;
					if (!chr.equals(flagchr)){
						LastPos=-1;
						startpoint=0;
						flagchr=chr;
					}
					for (int i=0; i< vec.size();i++){
						if (chr.equals(vec.get(i).ChrName)){
							for (int j= startpoint; j< vec.get(i).AlterStartVec.size();j++){
									if (position > vec.get(i).AlterEndVec.get(j) ){
										move=vec.get(i).AlterDisVec.get(j)  ;
										startpoint=j;
									}else{
										break;
									}
							}
						}
					}
					newpos+=move;
					if (newpos > LastPos ){
						LastPos = newpos;
						VcfWriter.write(chr+"\t");
						VcfWriter.write(newpos+"\t");
						for (int i=2;i<splitedInfor.length;i++ ){
							if (i!= (splitedInfor.length-1) ){
								VcfWriter.write(splitedInfor[i]+"\t");
							}else {
								VcfWriter.write(splitedInfor[i]+"\n");
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		VcfWriter.close();
	}
	void SomatinMuationOnGermlineInsertion(Vector <VcfChr> vec,String TumorVcfFile,  String TumorNewVcfFile, boolean removeduplicates) throws IOException{
		FileWriter mydata = new FileWriter(TumorNewVcfFile,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		FileReader fr =new FileReader(TumorVcfFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		int PrevPosition=-999; //(For Mutiple Contigs with Common Position, pgtool only reserve the first Position)
		//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
		VcfWriter.write("Chrom\t"+"Germ.Pos\t"+"Germ.Ref\t"+"Germ.Alt\t"+"Somatic.Rel.Pos\t"+"Somatic.Ref\t"+"Somatic.Alt\n");
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#"))  {
					String writenull=null;
					//VcfWriter.write(line+"\n")  ;   //Head Line
				}else {
					String[] splitedInfor = line.split("\t");
					String chr = splitedInfor[0].substring(0, splitedInfor[0].length());
					int position=Integer.parseInt(splitedInfor[1]);
					int newpos=-1;
					String tmp_ref=null;
					String tmp_alt=null;
					int tmp_pos=-1;
//					boolean haschr=false;
					for (int i=0; i< vec.size();i++){
						if (chr.equals(vec.get(i).ChrName)){
							int j=-1;
							int start= 0;
							int end=  vec.get(i).NewStartIndexVec.size()-1 ;
							while (start <= end) {
								int middle = (start + end) / 2;
								if ((position>= vec.get(i).NewStartIndexVec.get(middle))  && (position<= vec.get(i).NewEndIndexVec.get(middle))){
									j=middle;
									if (vec.get(i).InIndelIndexRegion.get(j) ){
										if ((!vec.get(i).InsertionRefVec.get(j).equals("")) &&  (!vec.get(i).InsertionAltVec.get(j).equals("")) ) {
											newpos= position-vec.get(i).NewStartIndexVec.get(j)+1;
											tmp_ref= vec.get(i).InsertionRefVec.get(j);
											tmp_alt= vec.get(i).InsertionAltVec.get(j);
											tmp_pos=vec.get(i).StartIndexVec.get(j);
											break;
										}else {
											newpos=-1;
											break;
										}

									}else{
										//newpos= position-vec.get(i).MovementIndexVec.get(j);
										newpos=-1;
										break;
									}
								}
								if (position < vec.get(i).NewStartIndexVec.get(middle) ) {
									 end = middle - 1;
								}else if (position >vec.get(i).NewEndIndexVec.get(middle)    ) {
						               start = middle + 1;
								}
							}
//							if (j==-1) {
//								System.out.println("error: can not map the somatic mutation to germline insertions\n");
//								System.out.println(line);
//							}

						}
					}
					if (newpos!=-1) {
						VcfWriter.write(chr+"\t");
						VcfWriter.write(tmp_pos+"\t");
						VcfWriter.write(tmp_ref+"\t");
						VcfWriter.write(tmp_alt+"\t");
						VcfWriter.write(newpos+"\t");
						if (splitedInfor.length>5) {
							VcfWriter.write(splitedInfor[3]+"\t");
							VcfWriter.write(splitedInfor[4]+"\t");
						}
						VcfWriter.write("\n");
					}

					if (((newpos!= PrevPosition ) && ( removeduplicates==true)) || ( removeduplicates==false)&& (1==0) ) {
						PrevPosition = newpos;
						VcfWriter.write(chr+"\t");
						VcfWriter.write(newpos+"\t");
						for (int i=2;i<splitedInfor.length;i++ ){
							if (i!= (splitedInfor.length-1) ){
								VcfWriter.write(splitedInfor[i]+"\t");
							}else {
								VcfWriter.write(splitedInfor[i]+"\n");
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		VcfWriter.close();
	}
	void MapVariants(Vector <VcfChr> vec,String TumorVcfFile,  String TumorNewVcfFile, boolean removeduplicates) throws IOException{
		FileWriter mydata = new FileWriter(TumorNewVcfFile,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		FileReader fr =new FileReader(TumorVcfFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		int PrevPosition=-999; //(For Mutiple Contigs with Common Position, pgtool only reserve the first Position)
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#"))  {
					VcfWriter.write(line+"\n")  ;   //Head Line
				}else {
					String[] splitedInfor = line.split("\t");
					String chr = splitedInfor[0].substring(0, splitedInfor[0].length());
					int position=Integer.parseInt(splitedInfor[1]);
					int newpos=position;
//					boolean haschr=false;
					for (int i=0; i< vec.size();i++){
						if (chr.equals(vec.get(i).ChrName)){
							int j=-1;
							int start= 0;
							int end=  vec.get(i).NewStartIndexVec.size()-1 ;
							while (start <= end) {
								int middle = (start + end) / 2;
								if ((position>= vec.get(i).NewStartIndexVec.get(middle))  && (position<= vec.get(i).NewEndIndexVec.get(middle))){
									j=middle;
									if (vec.get(i).InIndelIndexRegion.get(j) ){
										newpos=vec.get(i).MovementIndexVec.get(j);
										break;
									}else{
										newpos= position-vec.get(i).MovementIndexVec.get(j);
										break;
									}
								}
								if (position < vec.get(i).NewStartIndexVec.get(middle) ) {
									 end = middle - 1;
								}else if (position >vec.get(i).NewEndIndexVec.get(middle)    ) {
						               start = middle + 1;
								}
							}
//							if (j==-1) System.out.println("error: can not map the vcf file");
/*
							for (int j=0; j<vec.get(i).NewStartIndexVec.size();j++ ){
//								if ((chr.equals("chrM")) && (position>-1)){
//									System.out.println(position+" "+vec.get(i).NewStartIndexVec.size()+" "+vec.get(i).NewEndIndexVec.size()
//											+" "+vec.get(i).StartIndexVec.size()+" "+vec.get(i).EndIndexVec.size()+ vec.get(i).MovementIndexVec.size());
//									for (int s=0;s< vec.get(i).NewStartIndexVec.size();s++){
//										System.out.println(position+" "+vec.get(i).NewStartIndexVec.get(s)+" "+vec.get(i).NewEndIndexVec.get(s)
//												+" "+vec.get(i).StartIndexVec.get(s)+" "+vec.get(i).EndIndexVec.get(s)+" "+ vec.get(i).MovementIndexVec.get(s));
//									}
//								}
								if ((position>= vec.get(i).NewStartIndexVec.get(j))  && (position<= vec.get(i).NewEndIndexVec.get(j))){
									if (vec.get(i).InIndelIndexRegion.get(j) ){
										newpos=vec.get(i).MovementIndexVec.get(j);
										break;
									}else{
										newpos= position-vec.get(i).MovementIndexVec.get(j);
										break;
									}
								}
							}
*/
						}
					}
//					if (newpos!= -888 ){
					if (((newpos!= PrevPosition ) && ( removeduplicates==true)) || ( removeduplicates==false) ) {
						PrevPosition = newpos;
						VcfWriter.write(chr+"\t");
						VcfWriter.write(newpos+"\t");
						for (int i=2;i<splitedInfor.length;i++ ){
							if (i!= (splitedInfor.length-1) ){
								VcfWriter.write(splitedInfor[i]+"\t");
							}else {
								VcfWriter.write(splitedInfor[i]+"\n");
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		VcfWriter.close();
	}


	void AssignNewVcf(Vector <VcfChr> vec,String TumorVcfFile,  String TumorNewVcfFile) throws IOException{
		FileWriter mydata = new FileWriter(TumorNewVcfFile,true);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		FileReader fr =new FileReader(TumorVcfFile);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		int PrevPosition=-999; //(For Mutiple Contigs with Common Position, pgtool only reserve the first Position)
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("#"))  {
					VcfWriter.write(line+"\n")  ;   //Head Line
				}else {
					String[] splitedInfor = line.split("\t");
					String chr = splitedInfor[0].substring(0, splitedInfor[0].length());
					int position=Integer.parseInt(splitedInfor[1]);
					int newpos=position;
					for (int i=0; i< vec.size();i++){
						if (chr.equals(vec.get(i).ChrName)){
							int j=-1;
							int start= 0;
							int end=  vec.get(i).NewStartIndexVec.size()-1 ;
							while (start <= end) {
								int middle = (start + end) / 2;
								if ((position>= vec.get(i).NewStartIndexVec.get(middle))  && (position<= vec.get(i).NewEndIndexVec.get(middle))){
									j=middle;
									if (vec.get(i).InIndelIndexRegion.get(j) ){
										newpos=vec.get(i).MovementIndexVec.get(j);
										break;
									}else{
										newpos= position-vec.get(i).MovementIndexVec.get(j);
										break;
									}
								}
								if (position < vec.get(i).NewStartIndexVec.get(middle) ) {
									 end = middle - 1;
								}else if (position >vec.get(i).NewEndIndexVec.get(middle)    ) {
						               start = middle + 1;
								}
							}
							if (j==-1) System.out.println("error: can not map the vcf file");
/*
							for (int j=0; j<vec.get(i).NewStartIndexVec.size();j++ ){
//								if ((chr.equals("chrM")) && (position>-1)){
//									System.out.println(position+" "+vec.get(i).NewStartIndexVec.size()+" "+vec.get(i).NewEndIndexVec.size()
//											+" "+vec.get(i).StartIndexVec.size()+" "+vec.get(i).EndIndexVec.size()+ vec.get(i).MovementIndexVec.size());
//									for (int s=0;s< vec.get(i).NewStartIndexVec.size();s++){
//										System.out.println(position+" "+vec.get(i).NewStartIndexVec.get(s)+" "+vec.get(i).NewEndIndexVec.get(s)
//												+" "+vec.get(i).StartIndexVec.get(s)+" "+vec.get(i).EndIndexVec.get(s)+" "+ vec.get(i).MovementIndexVec.get(s));
//									}
//								}
								if ((position>= vec.get(i).NewStartIndexVec.get(j))  && (position<= vec.get(i).NewEndIndexVec.get(j))){
									if (vec.get(i).InIndelIndexRegion.get(j) ){
										newpos=vec.get(i).MovementIndexVec.get(j);
										break;
									}else{
										newpos= position-vec.get(i).MovementIndexVec.get(j);
										break;
									}
								}
							}
*/
						}
					}
					if (newpos!= -888 ){
//					if (newpos!= PrevPosition ){
						PrevPosition = newpos;
						VcfWriter.write(chr+"\t");
						VcfWriter.write(newpos+"\t");
						for (int i=2;i<splitedInfor.length;i++ ){
							if (i!= (splitedInfor.length-1) ){
								VcfWriter.write(splitedInfor[i]+"\t");
							}else {
								VcfWriter.write(splitedInfor[i]+"\n");
							}
						}
					}
				}
			}
		}
		bufferedreader.close();
		VcfWriter.close();
	}

	void RemoveOverlaps (Vector <String> ChrVec , Vector <VcfChr> v,  String output )throws IOException{
		FileWriter mydata = new FileWriter(output,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		VcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
		for (int i=0; i<v.size();i++ ){
			int pe=0;
			for (int j =0; j<v.get(i).StartPositionVec.size();j++ ){
				int start = v.get(i).StartPositionVec.get(j);
				int end= v.get(i).EndPositionVec.get(j);
//				System.out.println(start+"\t"+end +"\t"+pe);
				boolean isoverlap= false;
				if (start<=pe) isoverlap =true;
				if (!isoverlap){
					VcfWriter.write( v.get(i).ChrName+"\t"+ v.get(i).StartPositionVec.get(j) +"\t.\t"+
							v.get(i).RefStrVec.get(j)+"\t"+v.get(i).AltStrVec.get(j)+"\t.\t.\t.\tGT\t");
					if (v.get(i).IsHeteroVec.get(j)){
						VcfWriter.write("0/1\n");
					}else{
						VcfWriter.write("1/1\n");
					}
					if (end > pe){
						pe= end;
					}
				}else{
					if (end > pe){
						pe= end;
					}
				}
			}
		}
		VcfWriter.close();
	}

	void Select(Vector <VcfChr> v, String gt, String interval, String removeduplicates )throws IOException{

	}

	void SelectGenotype(String gt , Vector <VcfChr> v,  String output )throws IOException{
		FileWriter mydata = new FileWriter(output,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		VcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
		for (int i=0; i< v.size();i++){
			for (int j =0; j<v.get(i).StartPositionVec.size();j++ ){
				if (gt.equals("homo")){
					if (!v.get(i).IsHeteroVec.get(j)){
						VcfWriter.write( v.get(i).ChrName+"\t"+ v.get(i).StartPositionVec.get(j) +"\t.\t"+
								v.get(i).RefStrVec.get(j)+"\t"+v.get(i).AltStrVec.get(j)+"\t.\t.\t.\tGT\t") ;
						VcfWriter.write("1/1\n");
					}
				}else if (gt.equals("heter")){
					if (v.get(i).IsHeteroVec.get(j)){
						VcfWriter.write( v.get(i).ChrName+"\t"+ v.get(i).StartPositionVec.get(j) +"\t.\t"+
								v.get(i).RefStrVec.get(j)+"\t"+v.get(i).AltStrVec.get(j)+"\t.\t.\t.\tGT\t") ;
						VcfWriter.write("0/1\n");
					}
				}
			}
		}
		VcfWriter.close();
	}

	void SortVariants(Vector <String> ChrVec , Vector <VcfChr> v,  String output )throws IOException{
		FileWriter mydata = new FileWriter(output,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		VcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
		int totalnumofvars=0;
		for (int i=0; i< v.size();i++){
			totalnumofvars+= v.get(i).StartPositionVec.size();
		}
		for (int i=0; i< ChrVec.size();i++){
			int c=-1;
			String chr= ChrVec.get(i);
			for (int j=0; j< v.size();j++)
				if (v.get(j).ChrName.equals(chr) )
					c=j;
			if (c!=-1){
				for (int j =0; j<v.get(c).StartPositionVec.size();j++ ){
					for (int k =j; k<v.get(c).StartPositionVec.size();k++ ){
						if  (v.get(c).StartPositionVec.get(j)> v.get(c).StartPositionVec.get(k)){
							int tmp_int =-1;
							String tmp_str=null;
							boolean tmp_bool= false;
							tmp_int = v.get(c).StartPositionVec.get(j);
							v.get(c).StartPositionVec.set(j,  v.get(c).StartPositionVec.get(k) );
							v.get(c).StartPositionVec.set(k,tmp_int);
							tmp_str = v.get(c).RefStrVec.get(j);
							v.get(c).RefStrVec.set(j,  v.get(c).RefStrVec.get(k) );
							v.get(c).RefStrVec.set(k,tmp_str);
							tmp_str = v.get(c).AltStrVec.get(j);
							v.get(c).AltStrVec.set(j,  v.get(c).AltStrVec.get(k) );
							v.get(c).AltStrVec.set(k,tmp_str);
							tmp_bool = v.get(c).IsHeteroVec.get(j);
							v.get(c).IsHeteroVec.set(j,  v.get(c).IsHeteroVec.get(k) );
							v.get(c).IsHeteroVec.set(k,tmp_bool);
						}
					}
				}
				for (int j =0; j<v.get(c).StartPositionVec.size();j++ ){
					VcfWriter.write( v.get(c).ChrName+"\t"+ v.get(c).StartPositionVec.get(j) +"\t.\t"+
							v.get(c).RefStrVec.get(j)+"\t"+v.get(c).AltStrVec.get(j)+"\t.\t.\t.\tGT\t") ;
					if (v.get(c).IsHeteroVec.get(j)){
						VcfWriter.write("0/1\n");
					}else{
						VcfWriter.write("1/1\n");
					}
				}
			}
			System.out.println("Finish Sorting the Variants in Chromosome: "+chr);
		}
		VcfWriter.close();
	}

	void CombineVariants(Vector <String> ChrVec , Vector <VcfChr> v1, Vector <VcfChr> v2, String output )throws IOException{
		FileWriter mydata = new FileWriter(output,false);
		PrintWriter VcfWriter = new PrintWriter(mydata);
		VcfWriter.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
		for (int i=0; i< ChrVec.size();i++){
			int c1=-1, c2=-1;
			String chr= ChrVec.get(i);
			for (int j=0; j< v1.size();j++)
				if (v1.get(j).ChrName.equals(chr) )
					c1=j;
			for (int j=0; j< v2.size();j++)
				if (v2.get(j).ChrName.equals(chr) )
					c2=j;
			int p1=0, p2=0;
			if ( (c1!=-1 ) && (c2!=-1 )){
				while ((p1!= (v1.get(c1).StartPositionVec.size())  ) ||
						(p2!= (v2.get(c2).StartPositionVec.size())  )){
					if ((p1!=  (v1.get(c1).StartPositionVec.size())   )  &&
						(p2!=  (v2.get(c2).StartPositionVec.size()) ) ){
						if  (v1.get(c1).StartPositionVec.get(p1)<= v2.get(c2).StartPositionVec.get(p2) ){
							VcfWriter.write( v1.get(c1).ChrName+"\t"+ v1.get(c1).StartPositionVec.get(p1) +"\t.\t"+
									v1.get(c1).RefStrVec.get(p1)+"\t"+v1.get(c1).AltStrVec.get(p1)+"\t.\t.\t.\tGT\t") ;
							if (v1.get(c1).IsHeteroVec.get(p1)){
								VcfWriter.write("0/1\n");
							}else{
								VcfWriter.write("1/1\n");
							}
							p1++;
						}else{
							VcfWriter.write( v2.get(c2).ChrName+"\t"+ v2.get(c2).StartPositionVec.get(p2) +"\t.\t"+
									v2.get(c2).RefStrVec.get(p2)+"\t"+v2.get(c2).AltStrVec.get(p2)+"\t.\t.\t.\tGT\t") ;
							if (v2.get(c2).IsHeteroVec.get(p2)){
								VcfWriter.write("0/1\n");
							}else{
								VcfWriter.write("1/1\n");
							}
							p2++;
						}
					}else {
						if (p1!=  (v1.get(c1).StartPositionVec.size()) ){
							VcfWriter.write( v1.get(c1).ChrName+"\t"+ v1.get(c1).StartPositionVec.get(p1) +"\t.\t"+
									v1.get(c1).RefStrVec.get(p1)+"\t"+v1.get(c1).AltStrVec.get(p1)+"\t.\t.\t.\tGT\t") ;
							if (v1.get(c1).IsHeteroVec.get(p1)){
								VcfWriter.write("0/1\n");
							}else{
								VcfWriter.write("1/1\n");
							}
							p1++;
						}
						if (p2!=  (v2.get(c2).StartPositionVec.size()) ){
							VcfWriter.write( v2.get(c2).ChrName+"\t"+ v2.get(c2).StartPositionVec.get(p2) +"\t.\t"+
									v2.get(c2).RefStrVec.get(p2)+"\t"+v2.get(c2).AltStrVec.get(p2)+"\t.\t.\t.\tGT\t") ;
							if (v2.get(c2).IsHeteroVec.get(p2)){
								VcfWriter.write("0/1\n");
							}else{
								VcfWriter.write("1/1\n");
							}
							p2++;
						}
					}
				}
			}else{
				if (c1!=-1 ){
					for (p1=0; p1<v1.get(c1).StartPositionVec.size();p1++ ){
						VcfWriter.write( v1.get(c1).ChrName+"\t"+ v1.get(c1).StartPositionVec.get(p1) +"\t.\t"+
								v1.get(c1).RefStrVec.get(p1)+"\t"+v1.get(c1).AltStrVec.get(p1)+"\t.\t.\t.\tGT\t") ;
						if (v1.get(c1).IsHeteroVec.get(p1)){
							VcfWriter.write("0/1\n");
						}else{
							VcfWriter.write("1/1\n");
						}
					}
				}
				if (c2!=-1 ){
					for (p2=0; p2<v2.get(c2).StartPositionVec.size();p2++ ){
						VcfWriter.write( v2.get(c2).ChrName+"\t"+ v2.get(c2).StartPositionVec.get(p2) +"\t.\t"+
								v2.get(c2).RefStrVec.get(p2)+"\t"+v2.get(c2).AltStrVec.get(p2)+"\t.\t.\t.\tGT\t") ;
						if (v2.get(c2).IsHeteroVec.get(p2)){
							VcfWriter.write("0/1\n");
						}else{
							VcfWriter.write("1/1\n");
						}
						p2++;
					}
				}
			}
		}
		VcfWriter.close();
	}

	void Compare2Vcf(vcf v2) throws IOException{
		for (int i=0;i< VcfConfigVec.size();i++){
			boolean flag=false;
			for (int j=0; j< v2.VcfConfigVec.size();j++){
				if (v2.VcfConfigVec.get(j).ChrName.equals(VcfConfigVec.get(i).ChrName    )){
					if (v2.VcfConfigVec.get(j).StartPoint ==(VcfConfigVec.get(i).StartPoint    )){
//						System.out.println(v2.VcfConfigVec.get(j).ChrName+":"+ v2.VcfConfigVec.get(j).StartPoint );
					}
				}
			}
		}
	}

	void ReplaceSamGenotype(String samfile,String outputfile , int len) throws IOException{
		FileReader fr =new FileReader(samfile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		FileWriter mydata = new FileWriter(outputfile,false);
		PrintWriter fawrite = new PrintWriter(mydata);
		System.out.println( "Replacing Genetype \n");
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("@")){
					fawrite.write(line+"\n");
				}
				if (!line.substring(0, 1).equals("@")){
					String[] splited = line.split("\t");
					String chr= splited[2];
					if (chr.equals("*")){
						fawrite.write(line+"\n");
					}else{
						String chrq= chr.substring(0, chr.length() );
						int pos= Integer.parseInt(splited[3]);
						int start=-1, end=-1 , startindex=-1, endindex=-1, index=-1;
						boolean findchr=false;
						for (int i=0; i<VcfChrVec.size();i++ ){
							if (VcfChrVec.get(i).ChrName.equals( chrq)  ){
								findchr=true;
								start = 0;
								end = VcfChrVec.get(i).RefStrVec.size()-1;
								while (start <= end) {
									int middle = (start + end) / 2;
									if   (((VcfChrVec.get(i).StartPositionVec.get(middle) -pos )<  len  )  &&
									((VcfChrVec.get(i).EndPositionVec.get(middle) -pos)>-1 ) ) {
										index=middle;
										break;
									}
							        if (pos < VcfChrVec.get(i).StartPositionVec.get(middle) )  {
							               end = middle - 1;
							        } else if (pos > VcfChrVec.get(i).EndPositionVec.get(middle) )    {
							               start = middle + 1;
							        }
								}
								if (index!=-1){
									startindex=index;
									endindex=index;
									while (((VcfChrVec.get(i).StartPositionVec.get(startindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(startindex) -pos)>-1 ) ) {
										startindex--;
										if (startindex==-1 ) {
											break;
										}
									}
									while (((VcfChrVec.get(i).StartPositionVec.get(endindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(endindex) -pos)>-1 ) ){
										endindex++;
										if (endindex== VcfChrVec.get(i).EndPositionVec.size()){
											break;
										}
									}
									startindex++;
									endindex--;
									int[] startvec= new int[endindex-startindex+1] ;
									int[] endvec= new int[endindex-startindex+1] ;
									String[] refvec = new String[endindex-startindex+1] ;
									String[] altvec = new String[endindex-startindex+1] ;
									for (int j=0;j< startvec.length;j++){
										startvec[j]= VcfChrVec.get(i).StartPositionVec.get(startindex+j);
										endvec[j]= VcfChrVec.get(i).EndPositionVec.get(startindex+j);
										refvec[j] = VcfChrVec.get(i).RefStrVec.get(startindex+j);
										altvec[j]= VcfChrVec.get(i).AltStrVec.get(startindex+j);

									}
//									System.out.println( "The sam line is:  "+"\n"+line);
//									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
//											endvec, refvec ,altvec);
									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
											endvec, altvec, refvec);		// Remember that the alt and ref has changed

									for (int j=0; j<splited.length;j++ ){
										if ((j!=9) && (j!=10) && (j!= (splited.length-1 ) )){
											fawrite.write(splited[j] +"\t");
										}
										if (j==9){
											fawrite.write( tab9_10 +"\t");
										}
										if (j== (splited.length-1 ) ){
											fawrite.write(splited[j] +"\n");
										}
									}
								}else{
//									System.out.println("Can not find the corresponding refeneren genome:----------------------- "+line);
									fawrite.write(line+"\n");
								}
								break;
							}
						}
						if (!findchr) fawrite.write(line+"\n");
					}
				}
			}
		}
		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}

	void ReplaceSamGenotype_Meth(String samfile,String outputfile , int len) throws IOException{
		FileReader fr =new FileReader(samfile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		FileWriter mydata = new FileWriter(outputfile,false);
		PrintWriter fawrite = new PrintWriter(mydata);
		System.out.println( "Replacing Genetype \n");
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("@")){
					fawrite.write(line+"\n");
				}
				if (!line.substring(0, 1).equals("@")){
					String[] splited = line.split("\t");
					String chr= splited[2];
					if (chr.equals("*")){
						fawrite.write(line+"\n");
					}else{
						String chrq= chr.substring(0, chr.length() );
						int pos= Integer.parseInt(splited[3]);
						int start=-1, end=-1 , startindex=-1, endindex=-1, index=-1;
						boolean findchr=false;
//						for (int i=0; i<VcfChrVec.size();i++ ){
						if(chr_index_hm.containsKey(chrq)){
							int i=chr_index_hm.get(chrq);
							if (VcfChrVec.get(i).ChrName.equals( chrq)  ){
								findchr=true;
								start = 0;
								end = VcfChrVec.get(i).RefStrVec.size()-1;
								while (start <= end) {
									int middle = (start + end) / 2;
									if   (((VcfChrVec.get(i).StartPositionVec.get(middle) -pos )<  len  )  &&
									((VcfChrVec.get(i).EndPositionVec.get(middle) -pos)>-1 ) ) {
										index=middle;
										break;
									}
							        if (pos < VcfChrVec.get(i).StartPositionVec.get(middle) )  {
							               end = middle - 1;
							        } else if (pos > VcfChrVec.get(i).EndPositionVec.get(middle) )    {
							               start = middle + 1;
							        }
								}
								if (index!=-1){
									startindex=index;
									endindex=index;
									while (((VcfChrVec.get(i).StartPositionVec.get(startindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(startindex) -pos)>-1 ) ) {
										startindex--;
										if (startindex==-1 ) {
											break;
										}
									}
									while (((VcfChrVec.get(i).StartPositionVec.get(endindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(endindex) -pos)>-1 ) ){
										endindex++;
										if (endindex== VcfChrVec.get(i).EndPositionVec.size()){
											break;
										}
									}
									startindex++;
									endindex--;
									int[] startvec= new int[endindex-startindex+1] ;
									int[] endvec= new int[endindex-startindex+1] ;
									String[] refvec = new String[endindex-startindex+1] ;
									String[] altvec = new String[endindex-startindex+1] ;
									for (int j=0;j< startvec.length;j++){
										startvec[j]= VcfChrVec.get(i).StartPositionVec.get(startindex+j);
										endvec[j]= VcfChrVec.get(i).EndPositionVec.get(startindex+j);
										refvec[j] = VcfChrVec.get(i).RefStrVec.get(startindex+j);
										altvec[j]= VcfChrVec.get(i).AltStrVec.get(startindex+j);

									}
//									System.out.println( "The sam line is:  "+"\n"+line);
//									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
//											endvec, refvec ,altvec);
									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
											endvec, altvec, refvec);		// Remember that the alt and ref has changed

									for (int j=0; j<splited.length;j++ ){
										if ((j!=9) && (j!=10) && (j!= (splited.length-1 ) )){
											fawrite.write(splited[j] +"\t");
										}
										if (j==9){
											fawrite.write( tab9_10 +"\t");
										}
										if (j== (splited.length-1 ) ){
											fawrite.write(splited[j] +"\n");
										}
									}
								}else{
//									System.out.println("Can not find the corresponding refeneren genome:----------------------- "+line);
									fawrite.write(line+"\n");
								}
								break;
							}
						}
						if (!findchr) fawrite.write(line+"\n");
					}
				}
			}
		}
		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}

	void MultipleHits(String samfile,String outputfile ) throws IOException{
		Vector <String> sam_arr = new Vector<String>();
		Vector <String> f_chr_arr = new Vector<String>();
		Vector <String> f_pos_arr = new Vector<String>();
		Vector <String> s_chr_arr = new Vector<String>();
		Vector <String> s_pos_arr = new Vector<String>();

		FileReader fr =new FileReader(samfile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		FileWriter mydata = new FileWriter(outputfile,false);
		PrintWriter fawrite = new PrintWriter(mydata);
		String tmp_readname=null;
		boolean flag=false;
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("@")){
					fawrite.write(line+"\n");
				}else {
					String[] splited = line.split("\t");
					if (splited[0].equals(tmp_readname)){
						sam_arr.add(line);
						f_chr_arr.add(splited[2]);
						f_pos_arr.add(splited[3]);
						if (splited[6].equals("=")) {
							s_chr_arr.add(splited[2] );
						}else {
							s_chr_arr.add(splited[6] );
						}
						s_pos_arr.add(splited[7] );
					}else {
						flag=false;
						tmp_readname= splited[0];
						Vector <Integer> ok = new Vector<Integer>();
						ok.clear();
						int num= f_chr_arr.size();
						if (num==2) {
							flag=true;
							fawrite.write(sam_arr.get(0)+"\n");
							fawrite.write(sam_arr.get(1)+"\n");
						}else {
//							for (int i=0; i<num;i++ ) {
//								System.out.println(sam_arr.get(i) );
//								System.out.println(f_chr_arr.get(i) );
//								System.out.println(f_pos_arr.get(i) );
//								System.out.println(s_chr_arr.get(i) );
//								System.out.println(s_pos_arr.get(i) );
//							}
							for (int i=0; i<num;i++ ) {
								ok.add(0);
							}
							int index=0;
							for (int i=0; i<num;i++ ) {
								if (ok.get(i)==0) {
									for (int j=0; j<num;j++ ) {
										if ((f_chr_arr.get(i).equals( s_chr_arr.get(j)   ))  &&  (f_pos_arr.get(i).equals( s_pos_arr.get(j)   ))
											&& (s_chr_arr.get(i).equals( f_chr_arr.get(j)   ))   && (s_pos_arr.get(i).equals( f_pos_arr.get(j)   ))  && (i!=j)) {
											flag=true;
//											System.out.println(sam_arr.get(i) );
//											System.out.println(f_chr_arr.get(i) );
//											System.out.println(f_pos_arr.get(i) );
//											System.out.println(s_chr_arr.get(i) );
//											System.out.println(s_pos_arr.get(i) );
												index++;
												ok.set(i, 1);
												ok.set(j, 1);
												if (num!=1) {
													String[] ss = sam_arr.get(i).split("\t");
													ss[0]=ss[0]+":"+String.valueOf(index);
													for (int k=0; k<ss.length;k++ ) {
														if (k!= (ss.length-1) ) {
															fawrite.write(ss[k]+"\t");
														}else {
															fawrite.write(ss[k]+"\n");
														}
													}
													String[] tt = sam_arr.get(j).split("\t");
													tt[0]=tt[0]+":"+String.valueOf(index);
													for (int k=0; k<tt.length;k++ ) {
														if (k!= (tt.length-1) ) {
															fawrite.write(tt[k]+"\t");
														}else {
															fawrite.write(tt[k]+"\n");
														}
													}
												}
										}
									}
								}
							}
						}
//						if (flag==false) {
//							for (int i=0; i<sam_arr.size();i++ ) {
//								System.out.println(sam_arr.get(i) );
//							}
//						}
						sam_arr.clear();
				        f_pos_arr.clear();
				        f_chr_arr.clear();
				        s_pos_arr.clear();
				        s_chr_arr.clear();
				        sam_arr.add(line);
						f_chr_arr.add(splited[2]);
						f_pos_arr.add(splited[3]);
						if (splited[6].equals("=")) {
							s_chr_arr.add(splited[2] );
						}else {
							s_chr_arr.add(splited[6] );
						}
						s_pos_arr.add(splited[7] );

					}

				}
			}
		}

		bufferedreader.close();
		fawrite.flush();
		fawrite.close();

	}


	void MapSamRead(String samfile,String outputfile , int len) throws IOException{
		FileReader fr =new FileReader(samfile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line=null;
		FileWriter mydata = new FileWriter(outputfile,true);
		PrintWriter fawrite = new PrintWriter(mydata);
		System.out.println( "Replacing Genetype!\n");
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				if (line.substring(0, 1).equals("@")){
					fawrite.write(line+"\n");
				}
				if (!line.substring(0, 1).equals("@")){
					String[] splited = line.split("\t");
					String chr= splited[2];
					if (chr.equals("*")){
						fawrite.write(line+"\n");
					}else{
						String chrq= chr.substring(0, chr.length() );
						int pos= Integer.parseInt(splited[3]);
						int start=-1, end=-1 , startindex=-1, endindex=-1, index=-1;
						boolean findchr=false;
						for (int i=0; i<VcfChrVec.size();i++ ){
							if (VcfChrVec.get(i).ChrName.equals( chrq)  ){
								findchr=true;
								start = 0;
								end = VcfChrVec.get(i).RefStrVec.size()-1;
								while (start <= end) {
									int middle = (start + end) / 2;
									if   (((VcfChrVec.get(i).StartPositionVec.get(middle) -pos )<  len  )  &&
									((VcfChrVec.get(i).EndPositionVec.get(middle) -pos)>-1 ) ) {
										index=middle;
										break;
									}
							        if (pos < VcfChrVec.get(i).StartPositionVec.get(middle) )  {
							               end = middle - 1;
							        } else if (pos > VcfChrVec.get(i).EndPositionVec.get(middle) )    {
							               start = middle + 1;
							        }
								}
								if (index!=-1){
									startindex=index;
									endindex=index;
									while (((VcfChrVec.get(i).StartPositionVec.get(startindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(startindex) -pos)>-1 ) ) {
										startindex--;
										if (startindex==-1 ) {
											break;
										}
									}
									while (((VcfChrVec.get(i).StartPositionVec.get(endindex) -pos )<  len  )  &&
											((VcfChrVec.get(i).EndPositionVec.get(endindex) -pos)>-1 ) ){
										endindex++;
										if (endindex== VcfChrVec.get(i).EndPositionVec.size()){
											break;
										}
									}
									startindex++;
									endindex--;
									int[] startvec= new int[endindex-startindex+1] ;
									int[] endvec= new int[endindex-startindex+1] ;
									String[] refvec = new String[endindex-startindex+1] ;
									String[] altvec = new String[endindex-startindex+1] ;
									for (int j=0;j< startvec.length;j++){
										startvec[j]= VcfChrVec.get(i).StartPositionVec.get(startindex+j);
										endvec[j]= VcfChrVec.get(i).EndPositionVec.get(startindex+j);
										refvec[j] = VcfChrVec.get(i).RefStrVec.get(startindex+j);
										altvec[j]= VcfChrVec.get(i).AltStrVec.get(startindex+j);
									}

//									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
//											endvec, refvec ,altvec);

									String tab9_10= TransferRead( splited[9],splited[10],splited[5],pos,startvec,
											endvec, altvec, refvec);		// Remember that the alt and ref has changed
//									System.out.println( "The sam line is__________:  "+"\n"+line);
//									for (int j=0;j< startvec.length;j++){
//										System.out.println(startvec[j]+"\t"+ endvec[j]+"\t"+refvec[j]+"\t"+altvec[j]);
//									}
									for (int j=0; j<splited.length;j++ ){
										if ((j!=9) && (j!=10) && (j!= (splited.length-1 ) )){
											fawrite.write(splited[j] +"\t");
										}
										if (j==9){
											fawrite.write( tab9_10 +"\t");
										}
										if (j== (splited.length-1 ) ){
											fawrite.write(splited[j] +"\n");
										}
									}
								}else{
//									System.out.println("Can not find the corresponding refeneren genome:----------------------- "+line);
									fawrite.write(line+"\n");
								}
								break;
							}
						}
						if (!findchr) fawrite.write(line+"\n");
					}
				}
			}
		}
		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}
}
