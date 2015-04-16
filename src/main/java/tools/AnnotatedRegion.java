package tools;
import java.io.Serializable;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.StringUtils;

import util.Util;

public class AnnotatedRegion implements Serializable{
	public String annotation;
	public String chr;
	public int start;
	public int end;
	public int size;
	public char strand;
	public Map attributes;

	public AnnotatedRegion(String annotation, String chr, int start, int end, char strand) {
		this.chr = chr;
		this.annotation = annotation;
		this.start = start;
		this.end = end;
		this.strand = strand;
		size=end-start+1;
		attributes = new HashMap();
	}
	
	public AnnotatedRegion(String annotation, String chr, int start, int end, char strand, Map<String,Object> attributes) {
		this(annotation, chr, start, end, strand);
		this.attributes=attributes;
	}

	public int get5Prime(){
		Integer end5p = null;
		if(strand=='+')
			end5p = start;
		if(strand=='-')
			end5p = end;
		return end5p;
	}
	
	public int get5Prime(char strand){
		Integer end5p = null;
		if(strand=='+')
			end5p = start;
		if(strand=='-')
			end5p = end;
		return end5p;
	}
	
	public int get3Prime(){
		Integer end5p = null;
		if(strand=='-')
			end5p = start;
		if(strand=='+')
			end5p = end;
		return end5p;
	}
	
	public int get3Prime(char strand){
		Integer end3p = null;
		if(strand=='-')
			end3p = start;
		if(strand=='+')
			end3p = end;
		return end3p;
	}
	
	public int codingEnd(){
		return strand=='+' ? end : start;
	}
	
	public String toString(){
		return chr+":"+start+"-"+end;
	}
	
	public String toLongString(){
		return StringUtils.join(Util.list(annotation,strand,chr+":"+start+"-"+end),"\t");
	}
	
	public String toAttributeString(){
		List<String> attributesList = new LinkedList<String>();
		StringBuilder sb = new StringBuilder();
		for(Object key : attributes.keySet()){
			attributesList.add(key+"="+attributes.get(key));
		}
		return StringUtils.join(attributesList,";");
	}
	
	public String toGTFAttributeString(){
		List<String> attributesList = new LinkedList<String>();
		StringBuilder sb = new StringBuilder();
		for(Object key : attributes.keySet()){
			attributesList.add(key+" \""+attributes.get(key)+"\"");
		}
		return StringUtils.join(attributesList,";");
	}
	
	public static String GTFAttributeString(Map attributes){
		List<String> attributesList = new LinkedList<String>();
		
		for(Object key : attributes.keySet()){
			attributesList.add(key+" \""+attributes.get(key)+"\"");
		}
		
		return StringUtils.join(attributesList,";");
	}
	
	public static String attributeString(Map attributes, String key_value_format, String sep){
		List<String> attributesList = new LinkedList<String>();
		
		for(Object key : attributes.keySet()){
			attributesList.add(Util.sprintf(key_value_format, key, attributes.get(key)));
		}
		
		return StringUtils.join(attributesList,sep);
	}
	
	public static Map parseAttributeString(String attributes){
		Map attributeMap = new HashMap();
		for(String attribute : attributes.split(";")){
			String[] pair = attribute.split("=",2);
			attributeMap.put(pair[0], pair[1]);
		}
		return attributeMap;
	}
	
	public void addAttribute(String key, Object value){
		attributes.put(key, value);
	}
	
	public Object getAttribute(String key){
		return attributes.get(key);
	}

	public boolean isNegativeStrand() {
		return strand=='-';
	}
	
	public int hashCode() {
		return (toString()+strand).hashCode();
	}
	public boolean equals(Object o) {
		AnnotatedRegion other = (AnnotatedRegion) o;
		return other.start==start && other.end == end && other.strand==strand && other.chr.equals(chr);
	}
}