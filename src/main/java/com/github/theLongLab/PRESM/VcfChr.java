package com.github.theLongLab.PRESM;

import java.util.Vector;

public class VcfChr {
	public String ChrName;
	public Vector <Integer> StartPositionVec;
	public Vector <Integer> EndPositionVec;
	public Vector <Integer> AltEndPositionVec;
	public Vector <Integer> SourceVec;// Comes from Gatk or Pindel Vcf file
	public Vector <String> RefStrVec;
	public Vector <String> AltStrVec;
	public Vector <Boolean> IsHeteroVec;

	public Vector <Integer> AlterDisVec;
	public Vector <Integer> AlterStartVec;
	public Vector <Integer> AlterEndVec;

	public Vector <Integer> StartIndexVec;
	public Vector <Integer> EndIndexVec;
	public Vector <Integer> NewStartIndexVec ;
	public Vector <Integer> NewEndIndexVec;
	public Vector <Integer> MovementIndexVec;
	public Vector <Boolean> InIndelIndexRegion ;
	public Vector <String> InsertionRefVec;
	public Vector <String> InsertionAltVec;



	public boolean HasInterSection( int a, int b, int c, int d){
		if ((a>=c) && (a<=d)) return true;
		if ((b>=c) && (b<=d)) return true;
		if ((c>=a) && (c<=b)) return true;
		if ((d>=a) && (d<=b)) return true;
		return false;
	}


	public VcfChr(String chr, Vector <VcfContig> vec ){
		ChrName=chr;
		Vector <Integer> start_v = new Vector<Integer>();
		Vector <Integer> end_v = new Vector<Integer>();
		Vector <Integer> alt_end_v = new Vector<Integer>();
		Vector <String> ref_v = new Vector<String>();
		Vector <String> alt_v = new Vector<String>();
		Vector <Integer> flag_v = new Vector<Integer>();
		Vector <Boolean> gt_v = new Vector<Boolean>();

		for (int i=0; i<vec.size();i++ ){ //pindel vcf first{
			if (( vec.get(i).SourceFlag==1) && (vec.get(i).ChrName.equals(chr) )  ){
				start_v.add( vec.get(i).StartPoint );
				end_v.add(vec.get(i).EndPoint);
				alt_end_v.add(vec.get(i).AltEndPoint);
				ref_v.add(vec.get(i).RefStr);
				alt_v.add(vec.get(i).AltStr);
				flag_v.add(vec.get(i).SourceFlag );
				gt_v.add(vec.get(i).Ishetero );
			}
		}
		int pindelnum= start_v.size();

		for (int i=0; i<vec.size();i++ ){ //gatk vcf second
			if (( vec.get(i).SourceFlag==0) && (vec.get(i).ChrName.equals(chr) )  &&  (!vec.get(i).Filter.equals("lowQual") ) ){
				boolean in=false;
				for (int j=0;j< pindelnum  ;j++)
					if (HasInterSection(start_v.get(j),end_v.get(j), vec.get(i).StartPoint,vec.get(i).EndPoint ) ) in =true;
				if (!in){
					start_v.add( vec.get(i).StartPoint );
					end_v.add(vec.get(i).EndPoint);
					alt_end_v.add(vec.get(i).AltEndPoint);
					ref_v.add(vec.get(i).RefStr);
					alt_v.add(vec.get(i).AltStr);
					flag_v.add(vec.get(i).SourceFlag );
					gt_v.add(vec.get(i).Ishetero );

				}
			}
		}

		Vector <Integer> start_vv = new Vector<Integer>();
		Vector <Integer> end_vv = new Vector<Integer>();
		Vector <Integer> alt_end_vv = new Vector<Integer>();
		Vector <String> ref_vv = new Vector<String>();
		Vector <String> alt_vv = new Vector<String>();
		Vector <Integer> flag_vv = new Vector<Integer>();
		Vector <Boolean> gt_vv = new Vector<Boolean>();

		for (int i=0;i <start_v.size();i++ ){
			start_vv.add(start_v.get(i));
			end_vv.add(end_v.get(i));
			alt_end_vv.add(alt_end_v.get(i));
			ref_vv.add(ref_v.get(i));
			alt_vv.add(alt_v.get(i));
			flag_vv.add(flag_v.get(i));
			gt_vv.add(gt_v.get(i));
		}

		//Check the vector
		if (start_vv.size()>1){
			for (int i=1;i <start_vv.size();i++ ){
				if (HasInterSection( start_vv.get(i-1),end_vv.get(i-1),  start_vv.get(i),end_vv.get(i) ) ){
//				if (start_vv.get(i) <= end_vv.get(i-1)  ) {
					System.out.println( "Warning: There exists overlap between two variants: "+chr +":"
							+ start_vv.get(i)+"-"+end_vv.get(i) +" and "+chr +":"+  start_vv.get(i-1)+"-"+end_vv.get(i-1));
				}
			}
		}

		Vector <Integer> altervec = new Vector<Integer>();
		Vector <Integer> alterstartvec = new Vector<Integer>();
		Vector <Integer> alterendvec = new Vector<Integer>();

		altervec.add(0);
		alterstartvec.add(0);
		alterendvec.add(0);
		for (int i=0;i <start_vv.size();i++ ){
			if ( ref_vv.get(i).length()!= alt_vv.get(i).length() ){
				alterstartvec.add(start_vv.get(i)  );
				alterendvec.add(end_vv.get(i));
				int lengthofsize= altervec.size();
				altervec.add( altervec.get(lengthofsize-1) + alt_vv.get(i).length()- ref_vv.get(i).length()   );
			}
		}

		int lengthofsize= altervec.size();
		altervec.add( altervec.get(lengthofsize-1));
		alterstartvec.add(Integer.MAX_VALUE);
		alterendvec.add(Integer.MAX_VALUE);

		AlterDisVec =altervec ;
		AlterStartVec= alterstartvec;
		AlterEndVec= alterendvec;


//		System.out.println( chr+" indel:  "+AlterDisVec.get(AlterDisVec.size()-1 ));

		StartPositionVec= start_vv;
		EndPositionVec= end_vv;
		AltEndPositionVec= alt_end_vv;
		RefStrVec= ref_vv;
		AltStrVec= alt_vv;
		SourceVec= flag_vv;
		IsHeteroVec= gt_vv;

//		for (int i=0;i <StartPositionVec.size();i++ ){
//			System.out.println( StartPositionVec.get(i) +" "+ EndPositionVec.get(i)+" "+RefStrVec.get(i) );
//		}

// Construct Index for Vcf File
		Vector <Integer> StartVec = new Vector<Integer>();
		Vector <Integer> EndVec = new Vector<Integer>();
		Vector <Integer> NewStartVec = new Vector<Integer>();
		Vector <Integer> NewEndVec = new Vector<Integer>();
		Vector <Integer> MovementVec= new Vector<Integer>();
		Vector <Boolean> InIndelRegion = new Vector<Boolean>();
		Vector <Integer> RelativePosVec= new Vector<Integer>(); //0710
		Vector <String> AltInsVec= new Vector<String>(); //0710
		Vector <String> RefInsVec= new Vector<String>(); //0710
//Initialization

		int lastIndel=0;
		boolean hasindel=false;
		for (int i=0;i<start_vv.size();i++ ){
			if (ref_vv.get(i).length()!=alt_vv.get(i).length()  ){
				lastIndel=i;
				hasindel=true;
			}
		}
		if (hasindel){
			int move=0;
			StartVec.add(-1);
			NewStartVec.add(-1);
			InIndelRegion.add(false);
			AltInsVec.add(""); //0710
			RelativePosVec.add(0);//0710
			RefInsVec.add("");//0710

			MovementVec.add(0);
			for (int i=0;i<start_vv.size();i++ ){
				if (ref_vv.get(i).length()!=alt_vv.get(i).length()  ){

					if (i!=  lastIndel  ) {

						EndVec.add( start_vv.get(i)-1 );
						StartVec.add(start_vv.get(i));
						EndVec.add(end_vv.get(i));
						StartVec.add(end_vv.get(i)+1 );

						NewEndVec.add( start_vv.get(i)-1+move    );
						NewStartVec.add(start_vv.get(i)+move);
						move+=alt_vv.get(i).length()-ref_vv.get(i).length() ;
						NewEndVec.add(end_vv.get(i)+move);
						NewStartVec.add(end_vv.get(i)+1+move);

						MovementVec.add(start_vv.get(i)); //Retain the position of SNP
						MovementVec.add(move);

						InIndelRegion.add(true);
						InIndelRegion.add(false);
						if (alt_vv.get(i).length() > ref_vv.get(i).length() ) {
							AltInsVec.add(alt_vv.get(i));//0710
							RefInsVec.add(ref_vv.get(i));//0710
						}else {
							AltInsVec.add("");//0710
							RefInsVec.add("");//0710
						}
						AltInsVec.add("");//0710
						RefInsVec.add("");//0710

					}else{
						EndVec.add( start_vv.get(i)-1 );
						StartVec.add(start_vv.get(i));
						EndVec.add(end_vv.get(i));
						StartVec.add(end_vv.get(i)+1 );
						EndVec.add(1000000000 );

						NewEndVec.add( start_vv.get(i)-1+move    );
						NewStartVec.add(start_vv.get(i)+move);
						move+=alt_vv.get(i).length()-ref_vv.get(i).length() ;
						NewEndVec.add(end_vv.get(i)+move);
						NewStartVec.add(end_vv.get(i)+1+move);
						NewEndVec.add(1000000000+move );

						MovementVec.add(start_vv.get(i));
						MovementVec.add(move);

						InIndelRegion.add(true);
						InIndelRegion.add(false);

						if (alt_vv.get(i).length() > ref_vv.get(i).length() ) { //0710
							AltInsVec.add(alt_vv.get(i));//0710
							RefInsVec.add(ref_vv.get(i));//0710
						}else { //0710
							AltInsVec.add("");//0710
							RefInsVec.add("");//0710
						}
						AltInsVec.add("");//0710
						RefInsVec.add("");//0710
					}
				}

			}
		}

//		for (int i=0;i <StartVec.size();i++ ){
//			System.out.println("StartVec "+ StartVec.get(i)+" EndVec: "+EndVec.get(i)+" NewStartVec: "+NewStartVec.get(i)+" NewEndVec: "
//					+NewEndVec.get(i)+" MovementVec: "+MovementVec.get(i)+" InIndelRegion: "+InIndelRegion.get(i));
//		}

		StartIndexVec=StartVec;
		EndIndexVec=EndVec;
		NewStartIndexVec =NewStartVec;
		NewEndIndexVec=NewEndVec;
		MovementIndexVec=MovementVec;
		InIndelIndexRegion= InIndelRegion;
		InsertionRefVec=RefInsVec; //0710
		InsertionAltVec=AltInsVec; //0710
	}
}
