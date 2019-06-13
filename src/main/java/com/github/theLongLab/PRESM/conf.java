package com.github.theLongLab.PRESM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

public class conf {
	public Vector <String > ChrVec;
	public Vector <Integer > StartVec;
	public Vector <Integer > EndVec;
	public Vector <Integer > RefStartVec;
	public Vector <Integer > RefEndVec;
	public Vector <String > RefVec;
	public Vector <String > AltVec;

	public Vector <String> ChrNameVec;
	public Vector <Integer> ChrStartIndexVec;
	public Vector <Integer> ChrEndIndexVec;


	String transfer (String seq, int pos, int start, int end, int refstart, int refend, String ref, String alt, String flag) throws IOException{


		if (flag.equals("seq")) {
			System.out.println( seq+"\t"+ pos+"\t"+start+"\t"+end+"\t"+refstart+"\t"+refend+"\t"+ref+"\t"+alt);
		}
		refstart= refstart-pos;
		refend=refend-pos;
		if (refstart <0) refstart=0;
		if (refend <0) refend=-1;
		if ((refstart+1) > seq.length()) refstart = -1;
		if ((refend+1) > seq.length()) refend = seq.length()-1;
		StringBuilder sb = new StringBuilder(seq);

		if (flag.equals("seq")){
			if (( refend!=-1) && (refstart!=-1 )){

				int len = ref.length();
				if (len > sb.substring(refstart, refend+1).toString().length() ){
					len = sb.substring(refstart, refend+1).toString().length();
				}
				for (int i=0;i< len;i++){
					if  (! (  ref.substring(i, i+1).equals(  sb.substring(i, i+1).toString() ) ) ){
						return "####";
					}
				}
				System.out.println("real: "+sb.substring(refstart, refend+1).toString()+"\t"+ "ref: "+ ref);
				sb.replace(refstart, refend+1, alt);
			}
		}

		if (flag.equals("qual")){
			if (( refend!=-1) && (refstart!=-1 )){
				String qul= "";
				for (int i=0;i< alt.length();i++)
					qul=qul+"B";
				sb.replace(refstart, refend+1, qul);
			}
		}
		return sb.toString();
	}

	boolean noclip(String s1, String s2)throws IOException{
		int len =s1.length();
		for (int i=0; i< len;i++){
			if (s1.substring(i, i+1).equals(s2)  ){
				return false;
			}
		}
		return true;
	}



	void ReadConf(String file) throws IOException{

		Vector <String> tChrVec = new Vector<String>();
		Vector <String> tRefVec = new Vector<String>();
		Vector <String> tAltVec = new Vector<String>();
		Vector <Integer> tStartVec = new Vector<Integer>();
		Vector <Integer> tEndVec = new Vector<Integer>();
		Vector <Integer> tRefStartVec = new Vector<Integer>();
		Vector <Integer> tRefEndVec = new Vector<Integer>();

		FileReader fr =new FileReader(file);
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";
		while ( (line =bufferedreader.readLine())!=null ){
			if (0!=line.length()){
				String[] splited = line.split("\t");
				tChrVec.add(splited[0] );
				tStartVec.add(Integer.parseInt(splited[1]));
				tEndVec.add(Integer.parseInt(splited[2]));
				tRefStartVec.add(Integer.parseInt(splited[3]));
				tRefEndVec.add(Integer.parseInt(splited[4]));
				tRefVec.add(splited[5] );
				tAltVec.add(splited[6] );
			}
		}
		bufferedreader.close();
		Vector <String> tChrNameVec = new Vector<String>();
		Vector <Integer> tChrStartIndexVec = new Vector<Integer>();
		Vector <Integer> tChrEndIndexVec = new Vector<Integer>();

		tChrNameVec.add( tChrVec.get(0));
		tChrStartIndexVec.add( 0);
		String chr=tChrVec.get(0);
		for (int i=0; i<tChrVec.size();i++ ){
			if (!tChrVec.get(i).equals(chr)){
				tChrNameVec.add( tChrVec.get(i));
				tChrEndIndexVec.add(i-1);
				tChrStartIndexVec.add(i);
				chr= tChrVec.get(i);
			}
		}
		tChrEndIndexVec.add(tChrVec.size()-1);

//		for (int i=0;i<tChrVec.size();i++ )
//			System.out.println( tChrVec.get(i)+"\t"+ tStartVec.get(i)+"\t"+ tEndVec.get(i) + "\t"+ tRefStartVec.get(i) +
//					"\t"+ tRefEndVec.get(i) + "\t"+ tRefVec.get(i) + "\t"+ tAltVec.get(i));


		ChrNameVec= tChrNameVec;
		ChrStartIndexVec= tChrStartIndexVec;
		ChrEndIndexVec= tChrEndIndexVec;

		ChrVec= tChrVec;
		StartVec = tStartVec;
		EndVec= tEndVec;
		RefStartVec= tRefStartVec;
		RefEndVec= tRefEndVec;
		RefVec= tRefVec;
		AltVec= tAltVec;

	}

	void MapSamRead(String samfile,String outputfile ) throws IOException{
		FileReader fr =new FileReader(samfile );
		BufferedReader bufferedreader= new BufferedReader(fr);
		String line="";

		FileWriter mydata = new FileWriter(outputfile,true);
		PrintWriter fawrite = new PrintWriter(mydata);

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
						int start=0;
						int end=0;
						int index=-1;
						for (int i=0; i<ChrNameVec.size();i++ ){
							if (ChrNameVec.get(i).equals(chrq) ){
//								System.out.println(chr);
//								System.out.println (line);
//								System.out.println("position:  ");
								start= ChrStartIndexVec.get(i);
								end= ChrEndIndexVec.get(i);

							    while (start <= end) {
							    	int middle = (start + end) / 2;
							        if ((pos >= StartVec.get(middle)  ) && ( pos<=EndVec.get(middle)) ){
							        	index = middle;
							        	break;
							        }
							        if (pos < StartVec.get(middle)) {
							               end = middle - 1;
							        } else if (pos > EndVec.get(middle) ) {
							               start = middle + 1;
							        }
							    }
							}
						}

					    if (index==-1){
							   fawrite.write(line+"\n");
						}else{

							String tab9=transfer( splited[9],pos, StartVec.get( index),EndVec.get(index),
							    			RefStartVec.get(index), RefEndVec.get(index),RefVec.get(index), AltVec.get(index) ,"seq") ;
							String tab10=transfer( splited[10],pos, StartVec.get( index),EndVec.get(index),
							    			RefStartVec.get(index), RefEndVec.get(index),RefVec.get(index), AltVec.get(index) ,"qual") ;


							if ((noclip(splited[5],"S" ) )   && ( !tab9.equals("####")) ){
								for (int j=0; j<splited.length;j++ ){
									if ((j!=9) && (j!=10) && (j!= (splited.length-1 ) )){
										fawrite.write(splited[j] +"\t");
									}
									if (j==9){
										fawrite.write( tab9 +"\t");
									}
									if (j==10){
										fawrite.write( tab10 +"\t");
									}
									if (j== (splited.length-1 ) ){
										fawrite.write(splited[j] +"\n");
									}
								}
							}else {
//								System.out.println(line);
								fawrite.write(line+"\n");
							}

						}
					} //else if (!chr.substring(0, 1).equals("*"))

				}// if (!line.substring(0, 1).equals("@")){
			}// if (0!=line.length()){
		} //while

		bufferedreader.close();
		fawrite.flush();
		fawrite.close();
	}
}
