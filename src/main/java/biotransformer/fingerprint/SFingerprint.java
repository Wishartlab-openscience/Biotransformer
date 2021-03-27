/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.fingerprint;

import java.util.ArrayList;
import java.util.LinkedHashMap;

public class SFingerprint {
	
	private ArrayList<Double> bitValues 	= new  ArrayList<Double>();
	private double[] bitValuesDouble;
	private ArrayList<Integer> onbits	= new ArrayList<Integer>();
	private LinkedHashMap<String,Integer> featureMap = new LinkedHashMap<String,Integer>();

	public SFingerprint(ArrayList<Double> fingerprint_bits, LinkedHashMap<String, String> fpatterns) {
		// TODO Auto-generated constructor stub
		// TODO Auto-generated constructor stub
	this.bitValues = fingerprint_bits;
	this.collectOnBits();
	this.setBitValuesDouble();
	this.setFeatureMap(fpatterns);
	}


	public ArrayList<Double> getBitValues(){
		return this.bitValues;
	}
	public double[] getBitValuesDouble(){
		return this.bitValuesDouble;
	}
	public ArrayList<Integer> getOnBits(){
		return this.onbits;
	}

	public LinkedHashMap<String,Integer> getfeatureMap(){
		return this.featureMap;
	}
	
	
	public ArrayList<String> getOnFeatures(LinkedHashMap<String, String> fpatterns) throws Exception{
//		LinkedHashMap<String, String> fpatterns = ChemSearcher.getFingerprintPatterns();
		ArrayList<String> of = new ArrayList<String>();
		
		System.out.println(onbits.size());
		for(int j = 0; j < onbits.size(); j++){
			
			of.add((String) fpatterns.keySet().toArray()[onbits.get(j)-1]);
		}
		
		return of;
	}
	
	private void setFeatureMap(LinkedHashMap<String, String> fpatterns){
		String[] labels = fpatterns.keySet().toArray(new String[0]);
		for(int i = 0; i < this.getBitValues().size(); i++){
			this.featureMap.put(labels[i], this.getBitValues().get(i).intValue());
		}
	}

	private void collectOnBits(){
		for(int i =0; i < bitValues.size(); i++){

			if(bitValues.get(i)>=1.0){
				onbits.add(i+1);
			}
		}
	}

	private void setBitValuesDouble(){
		bitValuesDouble 	= new  double[bitValues.size()];
		for(int k =0; k < bitValues.size(); k++){
			bitValuesDouble[k] = bitValues.get(k);
			
		}
	}

}
