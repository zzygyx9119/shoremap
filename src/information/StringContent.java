package information;

import java.util.HashMap;

public class StringContent {

	private HashMap<String,HashMap<String, Integer>> counts= new HashMap<String,HashMap<String,Integer>>();
	
	
	public StringContent(){
		counts.put("stat", new HashMap<String, Integer>());
	}
	
	public boolean hasKey(String key){
		return !key.equals("stat") && counts.containsKey(key);
	}
	
	public HashMap<String, Integer> getCount(String key) throws Exception{
		if(key.equals("stat") || !counts.containsKey(key)){
			throw new Exception("Illegial key: "+key);
		}
		return counts.get(key);
	}
	
	public int getTotal(String key) throws Exception{
		if(key.equals("stat") || !counts.containsKey(key)){
			throw new Exception("Illegial key: "+key);
		}
		return counts.get("stat").get(key);
	}
	
	public void add(String a,String postfix){
		final int length= a.length();
		final String key=length+postfix;
		if(!counts.containsKey(key)){
			counts.put(key, new HashMap<String, Integer>());
			counts.get("stat").put(key, 0);
		}
		counts.get("stat").put(key, counts.get("stat").get(key)+1);
		if(!counts.get(key).containsKey(a)){
			counts.get(key).put(a, 1);
		}else{
			counts.get(key).put(a, counts.get(key).get(a)+1);
		}
	}
	
	public void add(String a){
		this.add(a, "");
	}
}
