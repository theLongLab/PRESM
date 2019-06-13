package com.github.theLongLab.PRESM;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

public class chrFasta{

	static final int LINE_LENGTH = 50;

	public Vector <String> ChrDnaVec = new Vector<String>();
	public Vector <Integer> ChrIndexVec = new Vector<Integer>();
	public String chrName;
	public String chrSeq;
	public List desc= new ArrayList();
	void assignment(String na, String seq) throws IOException{
		chrName=na.substring(0, na.length());
		chrSeq=seq;
	}

	boolean IsIUPACcode(String a){
		a=a.toUpperCase();
		if (a.equals("A")) return true;
		if (a.equals("N")) return true;
		if (a.equals("C")) return true;
		if (a.equals("G")) return true;
		if (a.equals("T")) return true;

		if (a.equals("U")) return true;
		if (a.equals("R")) return true;
		if (a.equals("Y")) return true;
		if (a.equals("S")) return true;
		if (a.equals("W")) return true;
		if (a.equals("K")) return true;
		if (a.equals("M")) return true;
		if (a.equals("B")) return true;
		if (a.equals("D")) return true;
		if (a.equals("H")) return true;
		if (a.equals("V")) return true;

//		if (a.equals(".")) return true;
//		if (a.equals("-")) return true;

		return false;
	}

	void MakePersonalizedReference(String file, Vector <VcfChr> vec, boolean flag) throws IOException{
		//First SNP and Deletion
				String seq= chrSeq;
//				System.out.println( seq.length()+"&&&&"  );
				seq="*"+seq; //* remember delete the space before writing the seq
				StringBuilder sb = new StringBuilder(seq);
				System.out.println( "Length of Chromosome:"+chrName +"  "+(sb.length()-1)  );
//				sb.replace(4, 7, "######");
				for (int i=0;i< vec.size();i++){
					if (vec.get(i).ChrName.equals(chrName)){
						for (int j=0;j< vec.get(i).AltStrVec.size();j++){
//							System.out.println( vec.get(i).AltStrVec.get(j).length() +"&&&"+ vec.get(i).RefStrVec.get(j).length());
//							System.out.println( vec.get(i).AltStrVec.get(j).length());
//							System.out.println( vec.get(i).RefStrVec.get(j).length());
						}
					}
				}


				Vector <String> SegmentSeq = new Vector<String>();

				Vector <String> HomoSegmentSeq = new Vector<String>();

				String HomoChr= "*"+chrName;
//				HomoSegmentSeq.add("NNNNNNNNNN");
//				String file2="ref.config";                  //0723 remove
//				FileWriter mydata2 = new FileWriter(file2,false);  //0723 remove
//				PrintWriter fawrite2 = new PrintWriter(mydata2);  //0723 remove



				for (int i=0;i< vec.size();i++){
//					System.out.println(  vec.get(i).ChrName  +" "+ chrName) ;
					if (vec.get(i).ChrName.equals(chrName)   ){

						Vector <Integer> startvec = new Vector<Integer>();
						Vector <Integer> endvec = new Vector<Integer>();
						Vector <String> refvec = new Vector<String>();
						Vector <String> altvec = new Vector<String>();
						Vector <Boolean> gtvec = new Vector<Boolean>();

						startvec= vec.get(i).StartPositionVec;
						endvec= vec.get(i).EndPositionVec;
						refvec=vec.get(i).RefStrVec;
						altvec =vec.get(i).AltStrVec;
						gtvec= vec.get(i).IsHeteroVec;

						Vector <Integer> gtstartvec = new Vector<Integer>();
						Vector <Integer> gtendvec = new Vector<Integer>();
						Vector <String> gtrefvec = new Vector<String>();
						Vector <String> gtaltvec = new Vector<String>();

						for (int j=0;j< startvec.size();j++){
							if (!gtvec.get(j)){
								gtstartvec.add(startvec.get(j));
								gtendvec.add(endvec.get(j)  );
								gtrefvec.add( refvec.get(j));
								gtaltvec.add(altvec.get(j)  );
							}
						}

//						long startindex= 11; //Remember Here is NOT 10 the first "NNNNNNNNNN" starts at pos 1 instead of 0

//						for (int j=0;j< startvec.size();j++){
//
//							if ((startvec.get(j)>= ReadLength ) && (endvec.get(j)<= ( seq.length()-ReadLength )  ) &&  (gtvec.get(j))   ){
//								if ((altvec.get(j).length()<ReadLength ) && (refvec.get(j).length()<ReadLength ) ){
//
//									fawrite2.write(chrName+"\t"+startindex+"\t"+(startindex +2*ReadLength-2+ refvec.get(j).length()-1 )+"\t"+
//									(startindex +ReadLength-1 )+"\t" +  (startindex +ReadLength-1 + refvec.get(j).length()-1)+"\t"+
//									refvec.get(j)+"\t"+altvec.get(j)+"\n") ;
//
//									String temp= seq.substring( startvec.get(j)- ReadLength+1, startvec.get(j) ) + refvec.get(j)+
//											seq.substring(endvec.get(j) +1,  endvec.get(j) + ReadLength )+"NNNNNNNNNN";
//
//									startindex =startindex+10+ 2*ReadLength-2+ refvec.get(j).length();
//
//									HomoSegmentSeq.add(temp);
//								}
//							}
//						}
						int len = seq.length();
						if ( gtstartvec.size()<1){
							HomoSegmentSeq.add( seq.substring(1, len) );
							break;
						}
						if (gtstartvec.get(0)!=1){
							HomoSegmentSeq.add( seq.substring(1,gtstartvec.get(0)  ));
							HomoSegmentSeq.add(gtaltvec.get(0));
						} else {

								HomoSegmentSeq.add( gtaltvec.get(0));
						}
						for (int j=1;j< gtstartvec.size();j++){
							HomoSegmentSeq.add(seq.substring(  (gtendvec.get(j-1)+1) ,  gtstartvec.get(j)   )  );
							HomoSegmentSeq.add(gtaltvec.get(j) );
						}
						if (gtendvec.get(gtendvec.size()-1)< seq.length())	HomoSegmentSeq.add( seq.substring(   (gtendvec.get(gtendvec.size()-1)+1),  seq.length() ));
					}


				}

//				fawrite2.flush();  //0723 remove
//				fawrite2.close();  //0723 remove

				boolean has=false;
				for (int i=0;i< vec.size();i++){

					if (vec.get(i).ChrName.equals(chrName)   ){
						has=true;
						Vector <Integer> startvec = new Vector<Integer>();
						Vector <Integer> endvec = new Vector<Integer>();
						Vector <String> refvec = new Vector<String>();
						Vector <String> altvec = new Vector<String>();

						startvec= vec.get(i).StartPositionVec;
						endvec= vec.get(i).EndPositionVec;
						refvec=vec.get(i).RefStrVec;
						altvec =vec.get(i).AltStrVec;
						int len = seq.length();

						if ( startvec.size()<1){
							SegmentSeq.add( seq.substring(1, len) );
							break;
						}
						if (startvec.get(0)!=1){
							SegmentSeq.add( seq.substring(1,startvec.get(0)  ));
							SegmentSeq.add(altvec.get(0));
						} else {
							SegmentSeq.add( altvec.get(0));
						}

						for (int j=1;j< startvec.size();j++){
//							if ((j % 10000)==0) System.out.println (chrName+" process "+i  );
							SegmentSeq.add(seq.substring(  (endvec.get(j-1)+1) ,  startvec.get(j)   )  );
							SegmentSeq.add(altvec.get(j) );
						}
						if (endvec.get(endvec.size()-1)< seq.length())	SegmentSeq.add( seq.substring(   (endvec.get(endvec.size()-1)+1),  seq.length() ));
					}
				}



				FileWriter mydata = new FileWriter(file,true);
				PrintWriter fawrite = new PrintWriter(mydata);
				fawrite.write(">"+chrName+"\n");
				for (int i=0;i< SegmentSeq.size();i++){
//					if ((i % 10000)==0) System.out.println (chrName+" finish "+i  );
					fawrite.write (SegmentSeq.get(i)  );
				}

				if (!has ) fawrite.write ( seq.substring(1, seq.length() )    );
//				if ((!flag) ||  ( HomoSegmentSeq.size()>0  ) ) fawrite.write("\n");
//
//				if ( HomoSegmentSeq.size()>0  ){
//					fawrite.write(">"+HomoChr+"\n");
//					for (int i=0;i< HomoSegmentSeq.size();i++){
//						fawrite.write (HomoSegmentSeq.get(i)  );
//					}
//				}

				if (!flag) fawrite.write("\n");



				fawrite.flush();
				fawrite.close();
		/*
				for (int i=0;i< vec.size();i++){
					if (vec.get(i).ChrName.equals(chrName)   ){
						for (int j=0;j< vec.get(i).AltStrVec.size();j++){
								if ((j % 10000 ) ==0)  System.out.println("snp "+chrName+" "+j+" "+ vec.get(i).StartPositionVec.get(j) );
								if ( vec.get(i).AltStrVec.get(j).length()==vec.get(i).RefStrVec.get(j).length() ){//SNP
									sb.replace(vec.get(i).StartPositionVec.get(j),  vec.get(i).EndPositionVec.get(j)+1,vec.get(i).AltStrVec.get(j));
								}
								if ( vec.get(i).AltStrVec.get(j).length()<vec.get(i).RefStrVec.get(j).length() ){ //deletion
									String temp= vec.get(i).AltStrVec.get(j);
									for (int s=0;s<(vec.get(i).RefStrVec.get(j).length()- vec.get(i).AltStrVec.get(j).length()) ;s++ ){
										temp=temp+"*";
									}

									sb.replace(vec.get(i).StartPositionVec.get(j),  vec.get(i).EndPositionVec.get(j)+1, temp);

							}
						}
					}
				}


		// Second Insertion

				for (int i=0;i< vec.size();i++){
					if (vec.get(i).ChrName.equals(chrName)   ){
						for (int j=(vec.get(i).AltStrVec.size()-1);j>-1;j--){
							if ((j % 10000 ) ==0)  System.out.println("insertion "+chrName+" "+j+" "+ vec.get(i).StartPositionVec.get(j) );
							if ( vec.get(i).AltStrVec.get(j).length()>vec.get(i).RefStrVec.get(j).length() ){ //Insrtion
								sb.replace( vec.get(i).StartPositionVec.get(j), vec.get(i).EndPositionVec.get(j)+1, vec.get(i).AltStrVec.get(j));
							}
						}
					}
				}
		//Last Step:Delete all stars in seq

				seq=sb.toString();
				seq=seq.replace("*","");
				System.out.println( "Length of New Chr:"+chrName +"  "+seq.length()  );


		// Write New Ref Fasta Format File
				FileWriter mydata = new FileWriter(file,true);
				PrintWriter fawrite = new PrintWriter(mydata);
				fawrite.write(">"+chrName+"\n");

				if (!flag) fawrite.write(seq+"\n");
				if (flag) fawrite.write(seq);


//				for (int i=120;(i<seq.length()-120);i++ ){
//					if ( !IsIUPACcode(seq.substring(i, i+1))  ) fawrite.write(seq.subSequence(i-100, i+100)+"  "+i+"\n");
//				}


//				long m=seq.length();
//				int nl= (int) (m /70);
//				for (int i=0;i< (nl+1);i++){
//					if (i!=nl){
//						fawrite.write(seq.substring(70*i, 70*i+70)+"\n");
//					}else{
//						if (!flag) fawrite.write(seq.substring(70*i,  (int)(m)  )+"\n");
//						if (flag ) fawrite.write(seq.substring(70*i,  (int)(m)  ) );
//					}
//				}

				fawrite.flush();
				fawrite.close();
		*/
			}


	void WriteNewFasta(String file, Vector <VcfChr> vec,  int ReadLength, boolean flag) throws IOException{
//First SNP and Deletion
		String seq= chrSeq;
//		System.out.println( seq.length()+"&&&&"  );
		seq="*"+seq; //* remember delete the space before writing the seq
		StringBuilder sb = new StringBuilder(seq);
		System.out.println( "Length of Hg38 Chr:"+chrName +"  "+(sb.length())  );
//		sb.replace(4, 7, "######");
		for (int i=0;i< vec.size();i++){
			if (vec.get(i).ChrName.equals(chrName)){
				for (int j=0;j< vec.get(i).AltStrVec.size();j++){
//					System.out.println( vec.get(i).AltStrVec.get(j).length() +"&&&"+ vec.get(i).RefStrVec.get(j).length());
//					System.out.println( vec.get(i).AltStrVec.get(j).length());
//					System.out.println( vec.get(i).RefStrVec.get(j).length());
				}
			}
		}


		Vector <String> SegmentSeq = new Vector<String>();

		Vector <String> HomoSegmentSeq = new Vector<String>();

		String HomoChr= "*"+chrName;
//		HomoSegmentSeq.add("NNNNNNNNNN");
		String file2="ref.config";
		FileWriter mydata2 = new FileWriter(file2,true);
		PrintWriter fawrite2 = new PrintWriter(mydata2);



		for (int i=0;i< vec.size();i++){
//			System.out.println(  vec.get(i).ChrName  +" "+ chrName) ;
			if (vec.get(i).ChrName.equals(chrName)   ){

				Vector <Integer> startvec = new Vector<Integer>();
				Vector <Integer> endvec = new Vector<Integer>();
				Vector <String> refvec = new Vector<String>();
				Vector <String> altvec = new Vector<String>();
				Vector <Boolean> gtvec = new Vector<Boolean>();

				startvec= vec.get(i).StartPositionVec;
				endvec= vec.get(i).EndPositionVec;
				refvec=vec.get(i).RefStrVec;
				altvec =vec.get(i).AltStrVec;
				gtvec= vec.get(i).IsHeteroVec;

				Vector <Integer> gtstartvec = new Vector<Integer>();
				Vector <Integer> gtendvec = new Vector<Integer>();
				Vector <String> gtrefvec = new Vector<String>();
				Vector <String> gtaltvec = new Vector<String>();

				for (int j=0;j< startvec.size();j++){
					if (!gtvec.get(j)){
						gtstartvec.add(startvec.get(j));
						gtendvec.add(endvec.get(j)  );
						gtrefvec.add( refvec.get(j));
						gtaltvec.add(altvec.get(j)  );
					}
				}

//				long startindex= 11; //Remember Here is NOT 10 the first "NNNNNNNNNN" starts at pos 1 instead of 0

//				for (int j=0;j< startvec.size();j++){
//
//					if ((startvec.get(j)>= ReadLength ) && (endvec.get(j)<= ( seq.length()-ReadLength )  ) &&  (gtvec.get(j))   ){
//						if ((altvec.get(j).length()<ReadLength ) && (refvec.get(j).length()<ReadLength ) ){
//
//							fawrite2.write(chrName+"\t"+startindex+"\t"+(startindex +2*ReadLength-2+ refvec.get(j).length()-1 )+"\t"+
//							(startindex +ReadLength-1 )+"\t" +  (startindex +ReadLength-1 + refvec.get(j).length()-1)+"\t"+
//							refvec.get(j)+"\t"+altvec.get(j)+"\n") ;
//
//							String temp= seq.substring( startvec.get(j)- ReadLength+1, startvec.get(j) ) + refvec.get(j)+
//									seq.substring(endvec.get(j) +1,  endvec.get(j) + ReadLength )+"NNNNNNNNNN";
//
//							startindex =startindex+10+ 2*ReadLength-2+ refvec.get(j).length();
//
//							HomoSegmentSeq.add(temp);
//						}
//					}
//				}
				int len = seq.length();
				if ( gtstartvec.size()<1){
					HomoSegmentSeq.add( seq.substring(1, len) );
					break;
				}
				if (gtstartvec.get(0)!=1){
					HomoSegmentSeq.add( seq.substring(1,gtstartvec.get(0)  ));
					HomoSegmentSeq.add(gtaltvec.get(0));
				} else {

						HomoSegmentSeq.add( gtaltvec.get(0));
				}
				for (int j=1;j< gtstartvec.size();j++){
					HomoSegmentSeq.add(seq.substring(  (gtendvec.get(j-1)+1) ,  gtstartvec.get(j)   )  );
					HomoSegmentSeq.add(gtaltvec.get(j) );
				}
				if (gtendvec.get(gtendvec.size()-1)< seq.length())	HomoSegmentSeq.add( seq.substring(   (gtendvec.get(gtendvec.size()-1)+1),  seq.length() ));
			}


		}

		fawrite2.flush();
		fawrite2.close();

		boolean has=false;
		for (int i=0;i< vec.size();i++){

			if (vec.get(i).ChrName.equals(chrName)   ){
				has=true;
				Vector <Integer> startvec = new Vector<Integer>();
				Vector <Integer> endvec = new Vector<Integer>();
				Vector <String> refvec = new Vector<String>();
				Vector <String> altvec = new Vector<String>();

				startvec= vec.get(i).StartPositionVec;
				endvec= vec.get(i).EndPositionVec;
				refvec=vec.get(i).RefStrVec;
				altvec =vec.get(i).AltStrVec;
				int len = seq.length();

				if ( startvec.size()<1){
					SegmentSeq.add( seq.substring(1, len) );
					break;
				}
				if (startvec.get(0)!=1){
					SegmentSeq.add( seq.substring(1,startvec.get(0)  ));
					SegmentSeq.add(altvec.get(0));
				} else {
					SegmentSeq.add( altvec.get(0));
				}

				for (int j=1;j< startvec.size();j++){
//					if ((j % 10000)==0) System.out.println (chrName+" process "+i  );
					SegmentSeq.add(seq.substring(  (endvec.get(j-1)+1) ,  startvec.get(j)   )  );
					SegmentSeq.add(altvec.get(j) );
				}
				if (endvec.get(endvec.size()-1)< seq.length())	SegmentSeq.add( seq.substring(   (endvec.get(endvec.size()-1)+1),  seq.length() ));
			}
		}



		FileWriter mydata = new FileWriter(file,true);
		PrintWriter fawrite = new PrintWriter(mydata);
		fawrite.write(">"+chrName+"\n");
		for (int i=0;i< SegmentSeq.size();i++){
//			if ((i % 10000)==0) System.out.println (chrName+" finish "+i  );
			fawrite.write (SegmentSeq.get(i)  );
		}

		if (!has ) fawrite.write ( seq.substring(1, seq.length() )    );
//		if ((!flag) ||  ( HomoSegmentSeq.size()>0  ) ) fawrite.write("\n");
//
//		if ( HomoSegmentSeq.size()>0  ){
//			fawrite.write(">"+HomoChr+"\n");
//			for (int i=0;i< HomoSegmentSeq.size();i++){
//				fawrite.write (HomoSegmentSeq.get(i)  );
//			}
//		}

		if (!flag) fawrite.write("\n");



		fawrite.flush();
		fawrite.close();
/*
		for (int i=0;i< vec.size();i++){
			if (vec.get(i).ChrName.equals(chrName)   ){
				for (int j=0;j< vec.get(i).AltStrVec.size();j++){
						if ((j % 10000 ) ==0)  System.out.println("snp "+chrName+" "+j+" "+ vec.get(i).StartPositionVec.get(j) );
						if ( vec.get(i).AltStrVec.get(j).length()==vec.get(i).RefStrVec.get(j).length() ){//SNP
							sb.replace(vec.get(i).StartPositionVec.get(j),  vec.get(i).EndPositionVec.get(j)+1,vec.get(i).AltStrVec.get(j));
						}
						if ( vec.get(i).AltStrVec.get(j).length()<vec.get(i).RefStrVec.get(j).length() ){ //deletion
							String temp= vec.get(i).AltStrVec.get(j);
							for (int s=0;s<(vec.get(i).RefStrVec.get(j).length()- vec.get(i).AltStrVec.get(j).length()) ;s++ ){
								temp=temp+"*";
							}

							sb.replace(vec.get(i).StartPositionVec.get(j),  vec.get(i).EndPositionVec.get(j)+1, temp);

					}
				}
			}
		}


// Second Insertion

		for (int i=0;i< vec.size();i++){
			if (vec.get(i).ChrName.equals(chrName)   ){
				for (int j=(vec.get(i).AltStrVec.size()-1);j>-1;j--){
					if ((j % 10000 ) ==0)  System.out.println("insertion "+chrName+" "+j+" "+ vec.get(i).StartPositionVec.get(j) );
					if ( vec.get(i).AltStrVec.get(j).length()>vec.get(i).RefStrVec.get(j).length() ){ //Insrtion
						sb.replace( vec.get(i).StartPositionVec.get(j), vec.get(i).EndPositionVec.get(j)+1, vec.get(i).AltStrVec.get(j));
					}
				}
			}
		}
//Last Step:Delete all stars in seq

		seq=sb.toString();
		seq=seq.replace("*","");
		System.out.println( "Length of New Chr:"+chrName +"  "+seq.length()  );


// Write New Ref Fasta Format File
		FileWriter mydata = new FileWriter(file,true);
		PrintWriter fawrite = new PrintWriter(mydata);
		fawrite.write(">"+chrName+"\n");

		if (!flag) fawrite.write(seq+"\n");
		if (flag) fawrite.write(seq);


//		for (int i=120;(i<seq.length()-120);i++ ){
//			if ( !IsIUPACcode(seq.substring(i, i+1))  ) fawrite.write(seq.subSequence(i-100, i+100)+"  "+i+"\n");
//		}


//		long m=seq.length();
//		int nl= (int) (m /70);
//		for (int i=0;i< (nl+1);i++){
//			if (i!=nl){
//				fawrite.write(seq.substring(70*i, 70*i+70)+"\n");
//			}else{
//				if (!flag) fawrite.write(seq.substring(70*i,  (int)(m)  )+"\n");
//				if (flag ) fawrite.write(seq.substring(70*i,  (int)(m)  ) );
//			}
//		}

		fawrite.flush();
		fawrite.close();
*/
	}
}


