package com.github.theLongLab.PRESM;

import java.lang.*;
import java.io.*;
import java.util.*;
import java.text.*;


public class CG {
	public static void main(String [] args) throws Exception {
		File dir = new File("");
		String CurrentPath= dir.getCanonicalPath().toString()+"/";
//		System.out.println( CurrentPath);
//      for(int i = 0; i < args.length; i++) {
//      	System.out.println(args[i]);
//      }
//	    String PindelVcfFile="/home/chencao/Workspace/data/test1.vcf";
//	    String PindelVcfFile="/home/chencao/Research/CancerGenome/RefGenome/NormalCompare2Hg38.v1/normal.pindel.vcf";
//	    String PindelVcfFile="/home/chencao/Research/CancerGenome/RefGenome/NormalCompare2Hg38.v1/normal.pindel.r2.vcf";
//	    String GatkVcfFile="/home/chencao/Research/CancerGenome/RefGenome/NormalCompare2Hg38.v1/normal.filter.vcf";
//	    String GatkVcfFile="/home/chencao/Workspace/data/test2.vcf";
//	    String VcfFile="/home/chencao/Workspace/data/ERR499.filter.vcf";
//	    String reffasta="/home/chencao/Research/CancerGenome/RefGenome/GRCh38.d1.vd1/GRCh38.d1.vd1.fa";
//	    String reffasta="/home/chencao/Research/CancerGenome/RefGenome/GRCh38.d1.vd1.new.fa";
//	    String reffasta="/home/chencao/Research/CancerGenome/RefGenome/test.fa";
//	    String OutputNewFa="/home/chencao/Research/CancerGenome/RefGenome/GRCh38.d1.vd1.new.fa";
//	    String BamFile="/home/chencao/Research/CancerGenome/RefGenome/tumor.bam";
//	    String OutputBamFile="/home/chencao/Research/CancerGenome/RefGenome/tumor.1.bam";
//	    Date date = new Date();
//      Human me= new Human();
//	    me.breath();

//		long cpustart =System.nanoTime();
//		long beforeUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//		System.out.println(beforeUsedMem);

		if(args.length<2){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n" +
			"Usage: java -jar presm.jar -F [function]");
			System.out.println("Supported functions:" +
			"\n\tCombineVariants" +
			"\n\tSortVariants" +
			"\n\tRemoveOverlaps" +
			"\n\tSelectGenotype" +
			"\n\tMakePersonalizedReference" +
			"\n\tMakePersonalizedVariantsDB" +
			"\n\tMapVariants"+
			"\n\tReplaceGenotype" +
			"\n\tSomaticMutationOnGermlineInsertion" +
			"\n\tViewFasta");
			System.exit(0);
		}

		String function=args[1];
		if(function.equals("ViewFasta")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
			"Developer: Chen Cao & Quan Long\n");
			System.out.println("View the specified region of the reference genome");
			if(args.length<3){
//				System.out.println("View the specified region of the reference genome");
				System.out.println("Usage: \n\t<-R\tinput reference genome file>\n\t" +
						"[-region\tspecified region]\n\t"+
						"[-L\tregion list]\n\t");
				System.exit(0);
			}else{
				String input=null;
				String region=null;
				String list= null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-R"))input=args[k+1];
						else if(args[k].equals("-region")) region=args[k+1];
						else if(args[k].equals("-L")) list=args[k+1];
					}
				}
				if(input==null){
					System.out.println("Input Fasta file can't be null!");
					System.exit(0);
				}
				if ((region==null) && (list==null) ) {
					System.out.println("Please specify the region or provide the region list file of the reference genome!");
					System.exit(0);
				}
				int min=1, max= Integer.MAX_VALUE;
				String chr=null;
				if (region!=null){
					String[] splited = region.split(":");
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
					fasta ff=new fasta( input);
					ff.ViewChr( chr,min,max);
					System.exit(0);
				}else if (list !=null ){
					FileReader fr =new FileReader(list);
					BufferedReader bufferedreader= new BufferedReader(fr);
					String line= null;
					fasta ff=new fasta( input);
					while ( (line =bufferedreader.readLine())!=null ){
						min=1;
						max= Integer.MAX_VALUE;
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

							ff.ViewChr( chr,min,max);
						}
					}
					bufferedreader.close();
					System.exit(0);
				}
			}
		}

		if(function.equals("CombineVariants")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Combine two variant records (VCF format) into one single vcf files.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-R\tinput reference genome file>\n\t" +
						"<-variant1\tinput1.vcf>\n\t"+
						"<-variant2\tinput2.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t");
			}else{
				String input1=null, input2=null, ref=null, output= null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-R"))ref=args[k+1];
						else if(args[k].equals("-variant1")) input1=args[k+1];
						else if(args[k].equals("-variant2")) input2=args[k+1];
						else if(args[k].equals("-O")) output=args[k+1];
					}
				}
				if(ref==null){
					System.out.println("Reference genome file can't be null!");
					System.exit(0);
				}
				if ((input1==null) || (input2==null) ) {
					System.out.println("Please designate the two variant files!");
					System.exit(0);
				}
				if(output==null){
					System.out.println("Please designate the output variant file!");
					System.exit(0);
				}
				fasta ff=new fasta( ref);
				vcf v1=new vcf();
				v1.readVcfFromFile(input1);
				vcf v2=new vcf();
				v2.readVcfFromFile(input2);
				vcf vc=new vcf();
				vc.CombineVariants(ff.ChrVec, v1.VcfChrVec, v2.VcfChrVec, output );
			}
			System.exit(0);
		}

		if(function.equals("SortVariants")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Sort variant file according to the reference genome.");
			if(args.length<3){

				System.out.println("Usage: \n\t<-R\tinput reference genome file>\n\t" +
						"<-variants\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t");
			}
			String input=null, ref=null, output=null;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-R"))ref=args[k+1];
					else if(args[k].equals("-variants")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
				}
			}
			if(ref==null){
				System.out.println("Reference genome file can't be null!");
				System.exit(0);
			}
			if (input==null ) {
				System.out.println("Please designate the variant file!");
				System.exit(0);
			}
			if(output==null){
				System.out.println("Please designate the output variant file!");
				System.exit(0);
			}
			fasta ff=new fasta( ref);
			vcf v=new vcf();
			v.readVcfFromFile(input);
			v.SortVariants(ff.ChrVec, v.VcfChrVec,output);
			System.exit(0);
		}

		if(function.equals("SelectGenotype")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
//			System.out.println("Select specified genotype of the variants.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-genotype\thomo/heter>\n\t" +
						"<-variants\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t");
			}
			String input=null, gt=null, output=null;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-genotype"))gt=args[k+1];
					else if(args[k].equals("-variants")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
				}
			}
			if(gt==null){
				System.out.println("Please designate the genotype!");
				System.exit(0);
			}
			if ((!gt.equals("homo")) && (!gt.equals("heter"))){
				System.out.println("Please designate the genotype: homo or heter!");
				System.exit(0);
			}
			if (input==null ) {
				System.out.println("Please designate the input variant file!");
				System.exit(0);
			}
			if(output==null){
				System.out.println("Please designate the output variant file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile(input);
			v.SelectGenotype(gt, v.VcfChrVec, output );
			System.exit(0);
		}

		if(function.equals("RemoveOverlaps")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Remove overlap variants in the vcf file.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-R\tinput reference genome file>\n\t" +
						"<-variants\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t");
			}
			String input=null, ref=null, output=null;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-R"))ref=args[k+1];
					else if(args[k].equals("-variants")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
				}
			}
			if(ref==null){
				System.out.println("Reference genome file can't be null!");
				System.exit(0);
			}
			if (input==null ) {
				System.out.println("Please designate the variant file!");
			}
			if(output==null){
				System.out.println("Please designate the output variant file!");
				System.exit(0);
			}
			fasta ff=new fasta( ref);
			vcf v=new vcf();
			v.readVcfFromFile(input);
			v.RemoveOverlaps(ff.ChrVec, v.VcfChrVec, output );
			System.exit(0);
		}

		if(function.equals("MakePersonalizedVariantsDB")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Make personalized variants database according to the germline mutaions provided by the user.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-variants\tvariant.vcf>\n\t" +
						"<-I\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tfile.intervals]\n\t"+
						"[-removeduplicates]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,variant=null;
			boolean removeduplicates=false;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-variants")) variant=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-removeduplicates")) removeduplicates=true;
				}
			}
			if(input==null){
				System.out.println("Please designate the input variant database file!");
				System.exit(0);
			}
			if (output==null ) {
				System.out.println("Please designate the output variant database file!");
				System.exit(0);
			}
			if(variant==null){
				System.out.println("Please designate the variant file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile(variant, gt, interval);

			v.MakePersonalizedVairants(v.VcfChrVec, input, output, removeduplicates);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
			System.exit(0);
		}

		if(function.equals("MakePersonalizedVariantsDB_Meth")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Make personalized variants database according to the germline mutaions provided by the user.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-variants\tvariant.vcf>\n\t" +
						"<-I\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tfile.intervals]\n\t"+
						"[-removeduplicates]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,variant=null;
			boolean removeduplicates=false;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-variants")) variant=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-removeduplicates")) removeduplicates=true;
				}
			}
			if(input==null){
				System.out.println("Please designate the input variant database file!");
				System.exit(0);
			}
			if (output==null ) {
				System.out.println("Please designate the output variant database file!");
				System.exit(0);
			}
			if(variant==null){
				System.out.println("Please designate the variant file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile_Meth(variant, gt, interval);
			v.MakePersonalizedVairants(v.VcfChrVec, input, output, removeduplicates);
			System.exit(0);
		}


		//--------------------------------------------------------------------------------------------
		//0710

		if(function.equals("SomaticMutationOnGermlineInsertion")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Output the relative coordinate of somatic mutations located on germline insertions.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-germlinemutations\tgermlinemutation.vcf>\n\t" +
						"<-I\tinput.vcf>\n\t"+
						"<-O\toutput.txt>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tinput.intervals]\n\t"+
						"[-removeduplicates]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,germlinemutation=null;
			boolean removeduplicates=false;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-germlinemutations"))germlinemutation=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-removeduplicates")) removeduplicates=true;
				}
			}
			if(input==null){
				System.out.println("Please designate the input variant file!");
				System.exit(0);
			}
			if (output==null ) {
				System.out.println("Please designate the output variant file!");
				System.exit(0);
			}
			if(germlinemutation==null){
				System.out.println("Please designate the germline mutation file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile(germlinemutation, gt, interval);
			v.SomatinMuationOnGermlineInsertion(v.VcfChrVec, input, output, removeduplicates);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
			System.exit(0);
		}

//0710
//---------------------------------------------------------------------------------------------
		if(function.equals("MapVariants")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Map variants to the universal reference genome according to the germline mutaions" +
					" provided by the users.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-germlinemutations\tgermlinemutation.vcf>\n\t" +
						"<-I\tinput.vcf>\n\t"+
						"<-O\toutput.vcf>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tinput.intervals]\n\t"+
						"[-removeduplicates]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,germlinemutation=null;
			boolean removeduplicates=false;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-germlinemutations"))germlinemutation=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-removeduplicates")) removeduplicates=true;
				}
			}
			if(input==null){
				System.out.println("Please designate the input variant file!");
				System.exit(0);
			}
			if (output==null ) {
				System.out.println("Please designate the output variant file!");
				System.exit(0);
			}
			if(germlinemutation==null){
				System.out.println("Please designate the germline mutation file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile(germlinemutation, gt, interval);
			v.MapVariants(v.VcfChrVec, input, output, removeduplicates);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
			System.exit(0);
		}

		if(function.equals("MakePersonalizedReference")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Build personalized reference genome according to the germline mutaions" +
					" provided by the user." );
			if(args.length<3){
				System.out.println("Usage: \n\t<-germlinemutations\tgermlinemutation.vcf>\n\t" +
						"<-I\tinput.fasta>\n\t"+
						"<-O\toutput.fasta>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tinput.intervals]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,germlinemutation=null;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-germlinemutations"))germlinemutation=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
				}
			}
			if(germlinemutation==null){
				System.out.println("Please designate the germline mutation file!");
				System.exit(0);
			}
			if(input==null){
				System.out.println("Please designate the input fasta file!");
				System.exit(0);
			}
			if (output==null ) {
				System.out.println("Please designate the output fasta file!");
				System.exit(0);
			}

			vcf v=new vcf();
			v.readVcfFromFile(germlinemutation, gt, interval);
			fasta ff=new fasta(input);
			for (int i=0;i< ff.chrFastaVec.size();i++){
		    	if (i!= (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).MakePersonalizedReference(output,v.VcfChrVec,false);
		    	if (i== (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).MakePersonalizedReference(output,v.VcfChrVec,true);
		    }
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
			System.exit(0);
		}

		if(function.equals("ReplaceGenotype")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Replace the alternative locus with reference locus accorindg to the germline heterzygous mutaions" +
					" provided by the users.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-germlinemutations\tgermline mutations.vcf>\n\t" +
						"<-I\tinput file>\n\t"+
						"<-O\toutput file>\n\t"+

						"<-readlength\tlength of read>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tfile.intervals]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,heterzygousmutations=null, filetype= null, ref=null;
			int readlength=-1;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-germlinemutations"))heterzygousmutations=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-readlength")) readlength =Integer.parseInt(args[k+1].toString());
				}
			}

			if(input==null){
				System.out.println("Please designate the input sequence alignment map file!");
				System.exit(0);
			}
			if (output==null){
				System.out.println("Please designate the output sequence alignment map file!");
				System.exit(0);
			}
			if(heterzygousmutations==null){
				System.out.println("Please designate the germline heterzygous mutation file!");
				System.exit(0);
			}
			if(readlength==-1){
				System.out.println("Please designate the length of the reads in the sequence alignment map file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile(heterzygousmutations, gt, interval, readlength);
			v.ReplaceSamGenotype(input, output,readlength) ;
//			long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//			long actualMemUsed=afterUsedMem-beforeUsedMem;
//			System.out.println(actualMemUsed);
//			long cputime = System.nanoTime() - cpustart;
//			System.out.printf("Each XXXXX took an average of %,d ns%n", cputime);
			System.exit(0);
		}

		if(function.equals("ReplaceGenotype_Meth")){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation discovery in cancer genomics\n" +
					"Developer: Chen Cao & Quan Long\n");
			System.out.println("Replace the alternative locus with reference locus accorindg to the germline heterzygous mutaions" +
					" provided by the users.");
			if(args.length<3){
				System.out.println("Usage: \n\t<-germlinemutations\tgermline mutations.vcf>\n\t" +
						"<-I\tinput file>\n\t"+
						"<-O\toutput file>\n\t"+

						"<-readlength\tlength of read>\n\t"+
						"[-genotype\thomo/heter]\n\t"+
						"[-intervals\tfile.intervals]\n\t");
			}
			String input=null, output=null, gt=null, interval=null,heterzygousmutations=null, filetype= null, ref=null;
			int readlength=-1;
			for(int k=1;k<args.length;k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-germlinemutations"))heterzygousmutations=args[k+1];
					else if(args[k].equals("-I")) input=args[k+1];
					else if(args[k].equals("-O")) output=args[k+1];
					else if(args[k].equals("-genotype")) gt=args[k+1];
					else if(args[k].equals("-intervals")) interval=args[k+1];
					else if(args[k].equals("-readlength")) readlength =Integer.parseInt(args[k+1].toString());
				}
			}

			if(input==null){
				System.out.println("Please designate the input sequence alignment map file!");
				System.exit(0);
			}
			if (output==null){
				System.out.println("Please designate the output sequence alignment map file!");
				System.exit(0);
			}
			if(heterzygousmutations==null){
				System.out.println("Please designate the germline heterzygous mutation file!");
				System.exit(0);
			}
			if(readlength==-1){
				System.out.println("Please designate the length of the reads in the sequence alignment map file!");
				System.exit(0);
			}
			vcf v=new vcf();
			v.readVcfFromFile_Meth(heterzygousmutations, gt, interval, readlength);
			v.ReplaceSamGenotype_Meth(input, output,readlength) ;
			System.exit(0);
		}


	    if (args[1].equals("MapRead")){
	    	String PindelVcfFile=args[3].toString();
	    	String GatkVcfFile= args[5].toString();
	    	String samfile= args[7].toString();
	    	String outputfile= args[9].toString();
	    	int ReadLength= Integer.parseInt(args[11].toString());
	    	ReadLength=ReadLength+20;
	    	vcf cc=new vcf();
	    	cc.readVcfFromFile(GatkVcfFile,null, null, ReadLength);
	    	cc.MapSamRead( samfile, outputfile,ReadLength) ;
	    	return;
	    }

	    if (args[1].equals("MultipleHits")){
	    	String SamFile=args[3].toString();
	    	String outputfile= args[5].toString();
	    	vcf cc=new vcf();
	    	cc.MultipleHits( SamFile, outputfile) ;
	    	return;
	    }


	    if (args[1].equals("MakeNewVcf")){
	    	String PindelVcfFile=args[3].toString();
	    	String GatkVcfFile= args[5].toString();
	    	vcf jj =new vcf();
	    	jj.readVcfFromFile(GatkVcfFile,PindelVcfFile);
	    	String TumorVcfFile= args[7].toString();
	    	String OutputNewVcfFile= args[9].toString();
	    	vcf jv =new vcf();
	    	jv.MakeNewVcf(jj.VcfChrVec ,TumorVcfFile,OutputNewVcfFile );
	    	return;
	    }

	    if (args[1].equals("ViewChr")){
	    	String reffasta=args[3].toString();
	    	String chr= args[5].toString();
	    	int start =Integer.parseInt(  args[7].toString());
	    	int end =Integer.parseInt(  args[9].toString());
	    	fasta ff=new fasta(reffasta );
	    	ff.ViewChr( chr,start,end);
	    }


	    if (args[1].equals("MakeFastq")){
	    	String reffasta=args[3].toString();
	    	String PindelVcfFile=args[5].toString();
	    	String GatkVcfFile= args[7].toString();
	    	vcf jj =new vcf();
	    	jj.readVcfFromFile(GatkVcfFile,PindelVcfFile);
	    	fasta ff=new fasta(reffasta );
	    	ff.MakeNewFastq(jj.VcfChrVec  );

	    }


	    if (false){
	    	String cigar="50M";
	    	String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
	    	int[] cigarnums = new int[cigarnumsTmp.length];
            for(int i = 0; i < cigarnumsTmp.length; i++){
            	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
            }
            String cigarlettersTmp[] = cigar.split("[0-9]+");//.toString().toCharArray();
            String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
        	for (int i=0;i< cigarletters.length;i++){
        		System.out.println(cigarletters[i] );
        		System.out.println(cigarnums[i] );
        	}
	    }

	    if (args[1].equals("MD")){
	    	String BamFile=args[3].toString() ;
	    	String OutputBamFile=args[5].toString();
		    Bam bb=new Bam();
		    bb.SamSelect(BamFile,OutputBamFile);
		    return ;
	    }

	    if (args[1].equals("FalseMapping")){
	    	String BamFile=args[3].toString() ;
	    	String OutputBamFile=args[5].toString();
		    Bam bb=new Bam();
		    bb.FalseMap(BamFile,OutputBamFile);
		    return ;
	    }

	    if (args[1].equals("SelectBamContig")){
	    	int SnpNumber=Integer.parseInt(  args[3].toString());
	    	String BamFile=args[5].toString() ;
	    	String OutputBamFile=args[7].toString();
		    Bam bb=new Bam();
		    bb.BamSelect(BamFile,OutputBamFile,SnpNumber);
		    return ;
	    }

	    if (args[1].equals("RevertSam")){
	    	String listfile= args[3].toString();
	    	String BamFile=args[5].toString() ;
	    	String OutputBamFile=args[7].toString();
		    Bam bb=new Bam();
		    bb.RevertSam(BamFile,OutputBamFile,listfile);
		    return ;
	    }

	    if (args[1].equals("RandomSplit")){
	    	int SplitNumber=Integer.parseInt(  args[5].toString());
	    	String Lib = args[3].toString();
	    	String BamFile=args[7].toString() ;
	    	String OutputBamFile=args[9].toString();
		    Bam bb=new Bam();
		    bb.SplitBam(BamFile,OutputBamFile,SplitNumber, Lib);
		    return ;
	    }

	    if (args[1].equals("RandomGenerateMutation")){
	    	int NumberOfMut=Integer.parseInt(  args[3].toString());
	    	String type = args[5].toString();
	    	String ref = args[9].toString();
	    	fasta ff=new fasta(ref);
	    	int len=Integer.parseInt(  args[7].toString());
	    	String OutputFile=args[11].toString();
	    	ff.RandomGenerateMuation(NumberOfMut, type, len, OutputFile);
	    }

	    if (args[1].equals("MakeNewFasta")){
	    	if (args.length==12){
		    	String reffasta=args[7].toString();
		    	String PindelVcfFile=args[3].toString();
		    	String GatkVcfFile= args[5].toString();
		    	String OutputNewFa=args[9].toString();
		    	int ReadLength= Integer.parseInt(args[11].toString());
		    	vcf jj =new vcf();
		    	jj.readVcfFromFile(GatkVcfFile,PindelVcfFile);
		    	fasta ff=new fasta(reffasta );
			    for (int i=0;i< ff.chrFastaVec.size();i++){
			    	if (i!= (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,ReadLength,false);
			    	if (i== (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,ReadLength,true );
			    }
	    	}
	    	if (args.length==10){
		    	String reffasta=args[5].toString();
		    	String GatkVcfFile= args[3].toString();
		    	String OutputNewFa=args[7].toString();
		    	int ReadLength= Integer.parseInt(args[9].toString());
		    	vcf jj =new vcf();
		    	jj.readVcfFromFile(GatkVcfFile);
		    	fasta ff=new fasta(reffasta );
			    for (int i=0;i< ff.chrFastaVec.size();i++){
			    	if (i!= (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,ReadLength,false);
			    	if (i== (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,ReadLength,true);
			    }
	    	}
	    	return;
	    }

	    if (args[1].equals("MapVcf")){
	    	String PindelVcfFile=args[3].toString();
	    	String GatkVcfFile= args[5].toString();
	    	vcf jj =new vcf();
	    	jj.readVcfFromFile(GatkVcfFile,PindelVcfFile);
	    	String TumorVcfFile= args[7].toString();
	    	String OutputNewVcfFile= args[9].toString();
//	    	int ReadLength= Integer.parseInt(args[9].toString());
	    	vcf jv =new vcf();
	    	jv.AssignNewVcf(jj.VcfChrVec ,TumorVcfFile,OutputNewVcfFile );
	    	return;
	    }
	    boolean printmessage=true;
	    if(printmessage){
			System.out.println("=============================================================================\n");
			System.out.println("PRESM1.0: Personalized Reference Editor for Somatic Mutation detection\n" +
			"Developer: Chen Cao & Quan Long\n" +
			"Usage: java -jar presm.jar -F [function]");
			System.out.println("Supported functions:" +
			"\n\tCombineVariants" +
			"\n\tSortVariants" +
			"\n\tRemoveOverlaps" +
			"\n\tSelectGenotype" +
			"\n\tMakePersonalizedReference" +
			"\n\tMakePersonalizedVariants" +
			"\n\tMapVariants"+
			"\n\tReplaceGenotype" +
			"\n\tViewFasta");
			System.exit(0);
		}

//	    if (args[1].equals("Simulation")){
//	    	if (args[3].equals("SNP")){
//	    		vcf jj =new vcf();
//	    		String GatkVcfFile = args[5].toString();
//	    		jj.readVcfFromFile(GatkVcfFile);
//
//	    		vcf db =new vcf();
//	    		String DbsnpVcfFile = args[7].toString();
//	    		db.readVcfFromFile(DbsnpVcfFile);
//	    	}
//	    }
//
//	    vcf jj =new vcf();
//	    jj.readVcfFromFile(GatkVcfFile,PindelVcfFile);
//	    System.out.println(jj.GatkVcfFilePath);
//
//	    fasta ff=new fasta(reffasta );
//
////	    System.out.println(ff.chrFastaVec.get(0).chrSeq.substring( 86087, 86096));
////		Sort the indel, snp, sv information in each chr
//
//
//	    int Mode=3;
//	    if (Mode==1){
//		    for (int i=0;i< ff.chrFastaVec.size();i++){
//		    	if (i!= (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,false);
//		    	if (i== (ff.chrFastaVec.size()-1)) ff.chrFastaVec.get(i).WriteNewFasta(OutputNewFa,jj.VcfChrVec,true);
//		    }
//	    }
//
//
//	    if (Mode==2){
//	    	String TumorVcfFile="/home/chencao/Research/CancerGenome/RefGenome/TumorCompare2Hg38.v1/tumor.merge.filter.r2.pass.vcf.recode.vcf";
//	    	String OutputNewVcfFile="/home/chencao/Desktop/tumor.merge.r2.new.pass.vcf";
//	    	vcf jv =new vcf();
//	    	jv.AssignNewVcf(jj.VcfChrVec ,TumorVcfFile,OutputNewVcfFile );
//	    }
//	    System.out.println("Game Over!!" );

//	    Hashtable<String,Integer> chr_index_hm  =new Hashtable<String,Integer>();
//		chr_index_hm.put("One", 111);
//		chr_index_hm.put("Two", 222);
//		chr_index_hm.put("Three", 333);
//		chr_index_hm.put("Four", 444);
//		if(chr_index_hm.containsKey("Two")){
//			System.out.println("Value: ");
//			System.out.println(chr_index_hm.get("Two"));
//		}
	}
}


class Human
{
    void breath()
    {
       System.out.println("hu...hu...ha...ha...die");
       height =181;
       age=18;
    }
    int height;
    public int age;
}
