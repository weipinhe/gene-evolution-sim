package com.webthinking.gene;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

public class Gene {
	private String name;

	private String sequence;
	private String orginalSequence;
	private String expressedProteinSeq;
	private String sequenceBeforeMutation;

	protected static Random randomGenerator = new Random(
			System.currentTimeMillis());

	public static boolean debugMode = false;

	public final static char BASAE_A = 'A', BASAE_T = 'T', BASAE_C = 'C',
			BASAE_G = 'G';
	public final static char[] BASE = { BASAE_A, BASAE_T, BASAE_C, BASAE_G };
	public final static int BASE_SIZE = BASE.length;

	public static Map<String, String> codonMap = new HashMap<String, String>(
			1000);

	static {
		// GCT, GCC, GCA, GCG to A
		codonMap.put("GCT", "A");
		codonMap.put("GCC", "A");
		codonMap.put("GCA", "A");
		codonMap.put("GCG", "A");

		// GAT, GAC to D
		codonMap.put("GAT", "D");
		codonMap.put("GAC", "D");

		// TTT, TTC to F
		codonMap.put("TTT", "F");
		codonMap.put("TTC", "F");

		// CAT, CAC to H
		codonMap.put("CAT", "H");
		codonMap.put("CAC", "H");

		// AAA, AAG
		codonMap.put("AAA", "K");
		codonMap.put("AAG", "K");

		// ATG
		codonMap.put("ATG", "M");

		// CCT, CCC, CCA, CCG
		codonMap.put("CCT", "P");
		codonMap.put("CCC", "P");
		codonMap.put("CCA", "P");
		codonMap.put("CCG", "P");

		// CGT, CGC, CGA, CGG, AGA, AGG
		codonMap.put("CGT", "R");
		codonMap.put("CGC", "R");
		codonMap.put("CGA", "R");
		codonMap.put("CGG", "R");
		codonMap.put("AGA", "R");
		codonMap.put("AGG", "R");

		// ACT, ACC, ACA, ACG
		codonMap.put("ACT", "T");
		codonMap.put("ACC", "T");
		codonMap.put("ACA", "T");
		codonMap.put("ACG", "T");

		// TGG
		codonMap.put("TGG", "W");

		// ATG
		// //!!!!!!!codonMap.put("ATG", "<");

		// TGT, TGC
		codonMap.put("TGT", "C");
		codonMap.put("TGC", "C");

		// GAA, GAG
		codonMap.put("GAA", "E");
		codonMap.put("GAG", "E");

		// GGT, GGC, GGA, GGG
		codonMap.put("GGT", "G");
		codonMap.put("GGC", "G");
		codonMap.put("GGA", "G");
		codonMap.put("GGG", "G");

		// ATT, ATC, ATA
		codonMap.put("ATT", "I");
		codonMap.put("ATC", "I");
		codonMap.put("ATA", "I");

		// TTA, TTG, CTT, CTC, CTA, CTG
		codonMap.put("TTA", "L");
		codonMap.put("TTG", "L");
		codonMap.put("CTT", "L");
		codonMap.put("CTC", "L");
		codonMap.put("CTA", "L");
		codonMap.put("CTG", "L");

		// AAT, AAC
		codonMap.put("AAT", "N");
		codonMap.put("AAC", "N");

		// CAA, CAG
		codonMap.put("CAA", "Q");
		codonMap.put("CAG", "Q");

		// TCT, TCC, TCA, TCG, AGT, AGC
		codonMap.put("TCT", "S");
		codonMap.put("TCC", "S");
		codonMap.put("TCA", "S");
		codonMap.put("TCG", "S");
		codonMap.put("AGT", "S");
		codonMap.put("AGC", "S");

		// GTT, GTC, GTA, GTG
		codonMap.put("GTT", "V");
		codonMap.put("GTC", "V");
		codonMap.put("GTA", "V");
		codonMap.put("GTG", "V");

		// TAT, TAC
		codonMap.put("TAT", "Y");
		codonMap.put("TAC", "Y");

		// TAA, TGA, TAG
		codonMap.put("TAA", ">");
		codonMap.put("TGA", ">");
		codonMap.put("TAG", ">");

		if (debugMode)
			System.out.println(codonMap);

	}

	//public static Map<String, Integer> AminoAcidScoreMap = new HashMap<String, Integer>(
	//		20);

	//static {

		//AminoAcidScoreMap.put("A", new Integer(1));
		//AminoAcidScoreMap.put("V", new Integer(2));
		//AminoAcidScoreMap.put("V", new Integer(3));
		//AminoAcidScoreMap.put("M", new Integer(20));
		//AminoAcidScoreMap.put("F", new Integer(6));
		//AminoAcidScoreMap.put("W", new Integer(10));
		//AminoAcidScoreMap.put("Y", new Integer(10));
		//AminoAcidScoreMap.put("V", new Integer(3));

		//AminoAcidScoreMap.put("S", new Integer(12));
		//AminoAcidScoreMap.put("T", new Integer(10));
		//AminoAcidScoreMap.put("N", new Integer(15));
		//AminoAcidScoreMap.put("M", new Integer(15));

		//AminoAcidScoreMap.put("C", new Integer(12));
		//AminoAcidScoreMap.put("U", new Integer(12));
		//AminoAcidScoreMap.put("G", new Integer(10));
		//AminoAcidScoreMap.put("P", new Integer(15));

		//AminoAcidScoreMap.put("M", new Integer(15));
		//AminoAcidScoreMap.put("M", new Integer(15));
		//AminoAcidScoreMap.put("M", new Integer(15));
		//AminoAcidScoreMap.put("D", new Integer(12));
		//AminoAcidScoreMap.put("M", new Integer(12));
		// AminoAcidScoreMap.put("A", new Integer(1));
		// AminoAcidScoreMap.put("A", new Integer(1));

	//}

	public static Map<String, Integer> AminoAcidSimilarityMap = new HashMap<String, Integer>(
			2000);

	// could be loaded from a config file
	static {
		// C combos
		AminoAcidSimilarityMap.put("CC", new Integer(12));
		AminoAcidSimilarityMap.put("CS", new Integer(0));
		AminoAcidSimilarityMap.put("CT", new Integer(-2));
		AminoAcidSimilarityMap.put("CP", new Integer(-3));
		AminoAcidSimilarityMap.put("CA", new Integer(-2));
		AminoAcidSimilarityMap.put("CG", new Integer(-3));
		AminoAcidSimilarityMap.put("CN", new Integer(-4));
		AminoAcidSimilarityMap.put("CD", new Integer(-5));

		AminoAcidSimilarityMap.put("CE", new Integer(-5));
		AminoAcidSimilarityMap.put("CQ", new Integer(-5));
		AminoAcidSimilarityMap.put("CH", new Integer(-3));
		AminoAcidSimilarityMap.put("CR", new Integer(-4));

		AminoAcidSimilarityMap.put("CK", new Integer(-5));
		AminoAcidSimilarityMap.put("CM", new Integer(-5));
		AminoAcidSimilarityMap.put("CI", new Integer(-2));
		AminoAcidSimilarityMap.put("CL", new Integer(-8));

		AminoAcidSimilarityMap.put("CV", new Integer(-2));
		AminoAcidSimilarityMap.put("CF", new Integer(-4));
		AminoAcidSimilarityMap.put("CY", new Integer(0));
		AminoAcidSimilarityMap.put("CW", new Integer(-8));

		// S combos
		AminoAcidSimilarityMap.put("SS", new Integer(2));
		AminoAcidSimilarityMap.put("ST", new Integer(1));
		AminoAcidSimilarityMap.put("SP", new Integer(1));
		AminoAcidSimilarityMap.put("SA", new Integer(1));
		AminoAcidSimilarityMap.put("SG", new Integer(1));
		AminoAcidSimilarityMap.put("SN", new Integer(1));
		AminoAcidSimilarityMap.put("SD", new Integer(0));
		AminoAcidSimilarityMap.put("SE", new Integer(0));

		AminoAcidSimilarityMap.put("SQ", new Integer(-1));
		AminoAcidSimilarityMap.put("SH", new Integer(-1));
		AminoAcidSimilarityMap.put("SR", new Integer(0));
		AminoAcidSimilarityMap.put("SK", new Integer(0));

		AminoAcidSimilarityMap.put("SM", new Integer(-2));
		AminoAcidSimilarityMap.put("SI", new Integer(-1));
		AminoAcidSimilarityMap.put("SL", new Integer(-3));
		AminoAcidSimilarityMap.put("SV", new Integer(-1));

		AminoAcidSimilarityMap.put("SF", new Integer(-3));
		AminoAcidSimilarityMap.put("SY", new Integer(-3));
		AminoAcidSimilarityMap.put("SW", new Integer(-2));

		// T combos
		AminoAcidSimilarityMap.put("TT", new Integer(3));
		AminoAcidSimilarityMap.put("TP", new Integer(0));
		AminoAcidSimilarityMap.put("TA", new Integer(1));
		AminoAcidSimilarityMap.put("TG", new Integer(0));
		AminoAcidSimilarityMap.put("TN", new Integer(0));
		AminoAcidSimilarityMap.put("TD", new Integer(0));
		AminoAcidSimilarityMap.put("TE", new Integer(0));

		AminoAcidSimilarityMap.put("TQ", new Integer(-1));
		AminoAcidSimilarityMap.put("TH", new Integer(-1));
		AminoAcidSimilarityMap.put("TR", new Integer(-1));
		AminoAcidSimilarityMap.put("TK", new Integer(0));

		AminoAcidSimilarityMap.put("TM", new Integer(-1));
		AminoAcidSimilarityMap.put("TI", new Integer(0));
		AminoAcidSimilarityMap.put("TL", new Integer(-2));
		AminoAcidSimilarityMap.put("TV", new Integer(0));

		AminoAcidSimilarityMap.put("TF", new Integer(-3));
		AminoAcidSimilarityMap.put("TY", new Integer(-3));
		AminoAcidSimilarityMap.put("TW", new Integer(-5));

		// P combos
		AminoAcidSimilarityMap.put("PP", new Integer(6));
		AminoAcidSimilarityMap.put("PA", new Integer(1));
		AminoAcidSimilarityMap.put("PG", new Integer(-1));
		AminoAcidSimilarityMap.put("PN", new Integer(-1));
		AminoAcidSimilarityMap.put("PD", new Integer(-1));
		AminoAcidSimilarityMap.put("PE", new Integer(-1));

		AminoAcidSimilarityMap.put("PQ", new Integer(0));
		AminoAcidSimilarityMap.put("PH", new Integer(0));
		AminoAcidSimilarityMap.put("PR", new Integer(0));
		AminoAcidSimilarityMap.put("PK", new Integer(-1));

		AminoAcidSimilarityMap.put("PM", new Integer(-2));
		AminoAcidSimilarityMap.put("PI", new Integer(-2));
		AminoAcidSimilarityMap.put("PL", new Integer(-3));
		AminoAcidSimilarityMap.put("PV", new Integer(-1));

		AminoAcidSimilarityMap.put("PF", new Integer(-5));
		AminoAcidSimilarityMap.put("PY", new Integer(-5));
		AminoAcidSimilarityMap.put("PW", new Integer(-6));

		// A combos
		AminoAcidSimilarityMap.put("AA", new Integer(2));
		AminoAcidSimilarityMap.put("AG", new Integer(1));
		AminoAcidSimilarityMap.put("AN", new Integer(0));
		AminoAcidSimilarityMap.put("AD", new Integer(0));
		AminoAcidSimilarityMap.put("AE", new Integer(0));

		AminoAcidSimilarityMap.put("AQ", new Integer(0));
		AminoAcidSimilarityMap.put("AH", new Integer(-1));
		AminoAcidSimilarityMap.put("AR", new Integer(-2));
		AminoAcidSimilarityMap.put("AK", new Integer(-1));

		AminoAcidSimilarityMap.put("AM", new Integer(-1));
		AminoAcidSimilarityMap.put("AI", new Integer(-1));
		AminoAcidSimilarityMap.put("AL", new Integer(-2));
		AminoAcidSimilarityMap.put("AV", new Integer(-0));

		AminoAcidSimilarityMap.put("AF", new Integer(-4));
		AminoAcidSimilarityMap.put("AY", new Integer(-3));
		AminoAcidSimilarityMap.put("AW", new Integer(-6));

		// G combos
		AminoAcidSimilarityMap.put("GG", new Integer(5));
		AminoAcidSimilarityMap.put("GN", new Integer(0));
		AminoAcidSimilarityMap.put("GD", new Integer(1));
		AminoAcidSimilarityMap.put("GE", new Integer(0));

		AminoAcidSimilarityMap.put("GQ", new Integer(-1));
		AminoAcidSimilarityMap.put("GH", new Integer(-2));
		AminoAcidSimilarityMap.put("GR", new Integer(-3));
		AminoAcidSimilarityMap.put("GK", new Integer(-2));

		AminoAcidSimilarityMap.put("GM", new Integer(-3));
		AminoAcidSimilarityMap.put("GI", new Integer(-3));
		AminoAcidSimilarityMap.put("GL", new Integer(-4));
		AminoAcidSimilarityMap.put("GV", new Integer(-1));

		AminoAcidSimilarityMap.put("GF", new Integer(-5));
		AminoAcidSimilarityMap.put("GY", new Integer(-5));
		AminoAcidSimilarityMap.put("GW", new Integer(-7));

		// N combos
		AminoAcidSimilarityMap.put("NN", new Integer(2));
		AminoAcidSimilarityMap.put("ND", new Integer(2));
		AminoAcidSimilarityMap.put("NE", new Integer(1));

		AminoAcidSimilarityMap.put("NQ", new Integer(1));
		AminoAcidSimilarityMap.put("NH", new Integer(2));
		AminoAcidSimilarityMap.put("NR", new Integer(0));
		AminoAcidSimilarityMap.put("NK", new Integer(1));

		AminoAcidSimilarityMap.put("NM", new Integer(-2));
		AminoAcidSimilarityMap.put("NI", new Integer(-2));
		AminoAcidSimilarityMap.put("NL", new Integer(-3));
		AminoAcidSimilarityMap.put("NV", new Integer(-2));

		AminoAcidSimilarityMap.put("NF", new Integer(-4));
		AminoAcidSimilarityMap.put("NY", new Integer(-2));
		AminoAcidSimilarityMap.put("NW", new Integer(-4));

		// D combos
		AminoAcidSimilarityMap.put("DD", new Integer(4));
		AminoAcidSimilarityMap.put("DE", new Integer(3));

		AminoAcidSimilarityMap.put("DQ", new Integer(2));
		AminoAcidSimilarityMap.put("DH", new Integer(1));
		AminoAcidSimilarityMap.put("DR", new Integer(-1));
		AminoAcidSimilarityMap.put("DK", new Integer(0));

		AminoAcidSimilarityMap.put("DM", new Integer(-3));
		AminoAcidSimilarityMap.put("DI", new Integer(-2));
		AminoAcidSimilarityMap.put("DL", new Integer(-4));
		AminoAcidSimilarityMap.put("DV", new Integer(-2));

		AminoAcidSimilarityMap.put("DF", new Integer(-6));
		AminoAcidSimilarityMap.put("DY", new Integer(-4));
		AminoAcidSimilarityMap.put("DW", new Integer(-7));

		// E combos
		AminoAcidSimilarityMap.put("EE", new Integer(4));

		AminoAcidSimilarityMap.put("EQ", new Integer(2));
		AminoAcidSimilarityMap.put("EH", new Integer(1));
		AminoAcidSimilarityMap.put("ER", new Integer(-1));
		AminoAcidSimilarityMap.put("EK", new Integer(0));

		AminoAcidSimilarityMap.put("EM", new Integer(-2));
		AminoAcidSimilarityMap.put("EI", new Integer(-2));
		AminoAcidSimilarityMap.put("EL", new Integer(-3));
		AminoAcidSimilarityMap.put("EV", new Integer(-2));

		AminoAcidSimilarityMap.put("EF", new Integer(-5));
		AminoAcidSimilarityMap.put("EY", new Integer(-4));
		AminoAcidSimilarityMap.put("EW", new Integer(-7));

		// Q combos
		AminoAcidSimilarityMap.put("QQ", new Integer(4));
		AminoAcidSimilarityMap.put("QH", new Integer(3));
		AminoAcidSimilarityMap.put("QR", new Integer(1));
		AminoAcidSimilarityMap.put("QK", new Integer(1));

		AminoAcidSimilarityMap.put("QM", new Integer(-1));
		AminoAcidSimilarityMap.put("QI", new Integer(-2));
		AminoAcidSimilarityMap.put("QL", new Integer(-2));
		AminoAcidSimilarityMap.put("QV", new Integer(-2));

		AminoAcidSimilarityMap.put("QF", new Integer(-5));
		AminoAcidSimilarityMap.put("QY", new Integer(-4));
		AminoAcidSimilarityMap.put("QW", new Integer(-5));

		// H combos
		AminoAcidSimilarityMap.put("HH", new Integer(6));
		AminoAcidSimilarityMap.put("HR", new Integer(2));
		AminoAcidSimilarityMap.put("HK", new Integer(0));

		AminoAcidSimilarityMap.put("HM", new Integer(-2));
		AminoAcidSimilarityMap.put("HI", new Integer(-2));
		AminoAcidSimilarityMap.put("HL", new Integer(-2));
		AminoAcidSimilarityMap.put("HV", new Integer(-2));

		AminoAcidSimilarityMap.put("HF", new Integer(-2));
		AminoAcidSimilarityMap.put("HY", new Integer(0));
		AminoAcidSimilarityMap.put("HW", new Integer(-3));

		// R combos
		AminoAcidSimilarityMap.put("RR", new Integer(8));
		AminoAcidSimilarityMap.put("RK", new Integer(3));

		AminoAcidSimilarityMap.put("RM", new Integer(0));
		AminoAcidSimilarityMap.put("RI", new Integer(-2));
		AminoAcidSimilarityMap.put("RL", new Integer(-3));
		AminoAcidSimilarityMap.put("RV", new Integer(-2));

		AminoAcidSimilarityMap.put("RF", new Integer(-4));
		AminoAcidSimilarityMap.put("RY", new Integer(-4));
		AminoAcidSimilarityMap.put("RW", new Integer(-2));

		// K combos
		AminoAcidSimilarityMap.put("KK", new Integer(5));

		AminoAcidSimilarityMap.put("KM", new Integer(0));
		AminoAcidSimilarityMap.put("KI", new Integer(-2));
		AminoAcidSimilarityMap.put("KL", new Integer(-3));
		AminoAcidSimilarityMap.put("KV", new Integer(-2));

		AminoAcidSimilarityMap.put("KF", new Integer(-5));
		AminoAcidSimilarityMap.put("KY", new Integer(-4));
		AminoAcidSimilarityMap.put("KW", new Integer(-3));

		// M combos
		AminoAcidSimilarityMap.put("MM", new Integer(6));
		AminoAcidSimilarityMap.put("MI", new Integer(2));
		AminoAcidSimilarityMap.put("ML", new Integer(4));
		AminoAcidSimilarityMap.put("MV", new Integer(2));

		AminoAcidSimilarityMap.put("MF", new Integer(0));
		AminoAcidSimilarityMap.put("MY", new Integer(-2));
		AminoAcidSimilarityMap.put("MW", new Integer(-4));

		// I combos
		AminoAcidSimilarityMap.put("II", new Integer(5));
		AminoAcidSimilarityMap.put("IL", new Integer(2));
		AminoAcidSimilarityMap.put("IV", new Integer(4));

		AminoAcidSimilarityMap.put("IF", new Integer(1));
		AminoAcidSimilarityMap.put("IY", new Integer(-1));
		AminoAcidSimilarityMap.put("IW", new Integer(-5));

		// L combos
		AminoAcidSimilarityMap.put("LL", new Integer(8));
		AminoAcidSimilarityMap.put("LV", new Integer(2));

		AminoAcidSimilarityMap.put("LF", new Integer(2));
		AminoAcidSimilarityMap.put("LY", new Integer(-1));
		AminoAcidSimilarityMap.put("LW", new Integer(-2));

		// V  combos
		AminoAcidSimilarityMap.put("VV", new Integer(4));

		AminoAcidSimilarityMap.put("VF", new Integer(-1));
		AminoAcidSimilarityMap.put("VY", new Integer(-2));
		AminoAcidSimilarityMap.put("VW", new Integer(-6));

		// F combos

		AminoAcidSimilarityMap.put("FF", new Integer(9));
		AminoAcidSimilarityMap.put("FY", new Integer(7));
		AminoAcidSimilarityMap.put("FW", new Integer(0));

		// Y combos

		AminoAcidSimilarityMap.put("YY", new Integer(10));
		AminoAcidSimilarityMap.put("YW", new Integer(0));

		// W combos

		AminoAcidSimilarityMap.put("WW", new Integer(17));
	}

	// CharSequence startCondo = ;

	public Gene(int size) {
		this(getRandomSequence(size));
	}

	public Gene(String inputSeq) {
		sequence = inputSeq;
		orginalSequence = sequence;
		if (debugMode) {
			this.outputSequence();
		}
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	protected static String getRandomSequence(int size) {
		StringBuffer randomSequence = new StringBuffer();
		int[] baseCount = new int[BASE_SIZE];
		/*
		 * for (int i = 0; i < size; i++) { int baseIndex =
		 * randomGenerator.nextInt(BASE_SIZE); if (debugMode)
		 * System.out.print(baseIndex); randomSequence.append(BASE[baseIndex]);
		 * baseCount[baseIndex]++;
		 *
		 * }
		 */

		int startIndex = 0;
		int baseIndex = 0;

		while (startIndex < size) {
			// endIndex = Math.min(startIndex + 3,sequence.length());
			String nextCodon = "";

			baseIndex = randomGenerator.nextInt(BASE_SIZE);
			nextCodon = nextCodon + BASE[baseIndex];
			baseCount[baseIndex]++;

			baseIndex = randomGenerator.nextInt(BASE_SIZE);
			nextCodon = nextCodon + BASE[baseIndex];
			baseCount[baseIndex]++;

			baseIndex = randomGenerator.nextInt(BASE_SIZE);
			nextCodon = nextCodon + BASE[baseIndex];
			baseCount[baseIndex]++;

			String AminoAcid = codonMap.get(nextCodon);
			if ("<".equals(AminoAcid) || ">".equals(AminoAcid))
				continue;
			else {
				randomSequence.append(nextCodon);
				startIndex = startIndex + 3;
			}
		}

		if (debugMode) {
			System.out.print('\n');
			System.out.print(size);
			System.out.print(';');
			for (int j = 0; j < baseCount.length; j++) {
				System.out.print((int) (baseCount[j] * 100 / size));
				System.out.print(';');
			}
			System.out.print('\n');

		}
		return randomSequence.toString();
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public String getOrginalSequence() {
		return orginalSequence;
	}

	public void setOrginalSequence(String orginalSequence) {
		this.orginalSequence = orginalSequence;
	}

	public void outputSequence() {
		System.out.println("Gene Name: " + name);
		System.out.println("DNA: " + sequence);
		System.out.println("PTN: " + express());
	}

	public void mutate(int times) {
		sequenceBeforeMutation = sequence;
		char[] newSequence = sequence.toCharArray();
		for (int i = 0; i < times; i++) {
			int position = randomGenerator.nextInt(sequence.length());
			newSequence[position] = BASE[randomGenerator.nextInt(BASE_SIZE)];
		}
		sequence = new String(newSequence);
		if (debugMode) {
			this.outputSequence();
		}
	}

	public void rejectMutate() {
		sequence = sequenceBeforeMutation;
	}

	public int compete(String sample) {
		int fitnessScore = 0;
		String proteinSeq = express();
		if (debugMode)
			System.out.println("------" + proteinSeq);
		if (proteinSeq.contains("<") || proteinSeq.contains(">")) {
			fitnessScore = -1000000;
			return fitnessScore;
		}

		fitnessScore = calculateFitnessScore(proteinSeq, sample);
		if (debugMode)
			System.out.println("------score = " + fitnessScore);
		return fitnessScore;
	}

	public int calculateFitnessScore(String proteinSeq, String sample) {
		int size = Math.min(proteinSeq.length(), sample.length());
		int fitnessScore = 0;
		for (int i = 0; i < size; i++) {
			String aminoAcidTarget = sample.substring(i, i + 1);
			String aminoAcidtoCompare = proteinSeq.substring(i, i + 1);


			Integer  positionSimilarityScore = AminoAcidSimilarityMap.get(aminoAcidTarget + aminoAcidtoCompare);
			if (positionSimilarityScore == null) {
				positionSimilarityScore = AminoAcidSimilarityMap.get(aminoAcidtoCompare + aminoAcidTarget );
			}
			int positionScore = positionSimilarityScore.intValue();

			fitnessScore += positionScore;

			//if (proteinSeq.charAt(i) == sample.charAt(i)) {
				//fitnessScore += positionScore;
			//} else {
				//String aminoAcidBad = proteinSeq.substring(i, i + 1);
				//int positionScoreBad = AminoAcidScoreMap.get(aminoAcidBad)
						//.intValue();
				//fitnessScore -= positionScoreBad;
				//;
			//}

		}
		//
		// int score = proteinSeq.hashCode();
		return fitnessScore;
	}

	public String express() {

		return transcribe(sequence);
	}

	public String transcribe(String dnaSequence) {
		StringBuffer proteinSequence = new StringBuffer();
		int startIndex = 0;
		int endIndex = 2;
		String nextCodon = "";
		while (startIndex < dnaSequence.length()) {
			endIndex = Math.min(startIndex + 3, dnaSequence.length());
			nextCodon = dnaSequence.substring(startIndex, endIndex);
			proteinSequence.append(codonMap.get(nextCodon));
			startIndex = startIndex + 3;
		}
		expressedProteinSeq = proteinSequence.toString();
		return expressedProteinSeq;
	}

	public int compareSequence(String seqToComapre) {
		return compareSequence(sequence, seqToComapre);
	}

	public int compareSequence(String sequence1, String sequence2) {
		char[] sequenceMark = sequence1.toCharArray();
		int changedCount = 0;
		for (int i = 0; i < sequence1.length(); i++) {
			sequenceMark[i] = (sequence1.charAt(i) == sequence2.charAt(i)) ? '='
					: '^';
			if (sequence1.charAt(i) != sequence2.charAt(i))
				changedCount++;
		}
		if (debugMode)
			System.out.println(sequence1);
		if (debugMode)
			System.out.println(sequence2);
		if (debugMode)
			System.out.println(new String(sequenceMark) + '!');
		if (debugMode)
			System.out.println("Difference: count = " + changedCount
					+ " out of " + sequence1.length() + "; rate = "
					+ (int) (changedCount * 100 / sequence1.length()) + "%");
		return changedCount;
	}

	public static Map<String, String> produceReverseProteinToDnaMap () {
		Map<String, String> proteinToDnaMap = new HashMap<String, String>();
		Set<Entry<String, String>> entrySet = codonMap.entrySet();
		Iterator<Entry<String, String>> itr = entrySet.iterator();
		while (itr.hasNext()) {
			String value = itr.next().getValue();
			String key = itr.next().getKey();
			proteinToDnaMap.put(value, key);
		};
		
		return proteinToDnaMap;
	}
	
	public static String decodeProteinToDNA (String proteiSequence) {
		Map<String, String> proteinToDnaMap = produceReverseProteinToDnaMap();
		
		StringBuffer dnaSequence = new StringBuffer();
		for (int i = 0; i < proteiSequence.length(); i++) {
			String currentAminoAcidCode  = proteiSequence.substring(i, i+1);
			dnaSequence.append(proteinToDnaMap.get(currentAminoAcidCode));
			
		}
		return dnaSequence.toString();
	}

	public int compareWithOriginalSequence() {
		return this.compareSequence(this.getOrginalSequence());
	}

	public static void main(String[] args) {
		Gene.debugMode = true;

		if (Gene.debugMode)
			System.out.println(Gene.codonMap);

		testGeneCreation();

		evolve();

	}

	public static void testGeneCreation() {
		// Compare two randomly generated genes' sequence
		System.out
				.println("\n******** Compare two randomly generated genes' sequence: ");
		Gene randomGene1 = new Gene(2000);
		randomGene1.compete("");
		Gene randomGene2 = new Gene(2000);
		randomGene2.compete("");
		// int randomChangedCount =
		// randomGene1.compareSequence(randomGene2.getSequence());
	}

	public static void evolve() {
		// ********* start testing mutation
		Gene.debugMode = true;
		System.out.println("\n");
		System.out.println("******** Start testing mutation and evolve:");
		System.out.println("Input needed!!! Copy and paste in your target protein sequence (type in \"skip\" if you want to use the default):");
		
		Scanner scaner = new Scanner(System.in);

		String targetProteinSeq = // "MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESFTFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLY"
									// + "";
		"FNVILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDPNEVFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPETSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYVYLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYGDVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQDHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSANATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIACGVAIPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLNNPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQSQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDSYGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRVDFSNKSNVNVGQVKDIHGRIPEML";

		String inputTargetSequence = scaner.next(); 
		if (!"skip".equalsIgnoreCase(inputTargetSequence)) {
			targetProteinSeq = inputTargetSequence;
		}
		
		int generation = 1000000;
		
		System.out.println("Generations to evolve (type 0 to select the default value of 1,000,000): ");		
		int inputGeneration = scaner.nextInt();
		if (inputGeneration > 0 ) {
			generation = inputGeneration;
		}
		
		Gene testGene;
		System.out.println("Sart sequence (type \"skip\" to select randomly generated sequnce of the same length): ");
		String inputStartSequence = scaner.next(); 
		if ("skip".equalsIgnoreCase(inputStartSequence)) {
			testGene = new Gene(targetProteinSeq.length() * 3);
		}
		else {
			String startDnaSeq = decodeProteinToDNA(inputStartSequence);
			testGene = new Gene(startDnaSeq);
		}
		
		scaner.close();
		
		// testGene.compete("");

		Gene.debugMode = false;
		long startTime = System.currentTimeMillis();

		int mutationStrength = 3;
		int selectionPressureTolerance = 1; // smaller is bigger pressure. 0
											// being biggest pressure
		int i, changedCount, score, oldScore = -1000000, rejectCount = 0, lastScoreChangeIndex = 0;
		for (i = 0; i < generation; i++) {
			testGene.mutate(mutationStrength);
			score = testGene.compete(targetProteinSeq);
			if (score < (oldScore - randomGenerator.nextInt(selectionPressureTolerance))
					|| score == -1000000) {
				testGene.rejectMutate();
				rejectCount++;
			} else {
				if (score != oldScore) {
					lastScoreChangeIndex = i;
				}
				oldScore = score;
			}
			int interval = 10000;
			if (i < interval) interval = 100; 
			if (i % interval == 0) {
				if (Gene.debugMode)
					System.out.println("++Mutation count = " + i);
				changedCount = testGene.compareWithOriginalSequence();
				System.out.println(i + "," + changedCount + "," + oldScore);
			}
		}

		long endTime = System.currentTimeMillis();
		long timeUsed = endTime - startTime;

		System.out.println("\nFinal mutation result: " + "gene size = "
				+ testGene.getSequence().length() + "; mutation strength = "
				+ mutationStrength + "; generation = " + i
				+ "; reject count = " + rejectCount + "; \n\tfinal score = "
				+ oldScore + "; last score change index = "
				+ lastScoreChangeIndex + "; Time used is " + timeUsed + " ms.");

		Gene.debugMode = true;
		// int finalChangedCount =
		testGene.compareWithOriginalSequence();
		testGene.compete("");

		System.out.println("\n**** Compare with orginal protein sequence:");
		testGene.compareSequence(
				testGene.transcribe(testGene.getOrginalSequence()),
				testGene.transcribe(testGene.getSequence()));

		System.out.println("**** Compare with target protein sequence:");
		testGene.compareSequence(targetProteinSeq,
				testGene.transcribe(testGene.getSequence()));

		// ********* end testing mutation
	}
}
