package tools.shoreMap;

import java.util.ArrayList;
import java.util.Collection;

public class ShoreMapLines {

	ArrayList<ShoreMapLine> lines;
	
	public ShoreMapLines(){
		lines=new ArrayList<ShoreMapLine>();
	}
	
	public ShoreMapLines(ShoreMapLine sml){
		this();
		lines.add(sml);
	}
	
	public ShoreMapLines(Collection<ShoreMapLine> l){
		this();
		lines.addAll(l);
	}
}
