package utils;

import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;

public class LineReader {

	private final static int bufferSize=2097152;//65536;
//	private StringBuffer strBuf;
	private FileChannel inChannel;
//	private ByteBuffer buf;
	
	private ByteBuffer buf= ByteBuffer.allocateDirect(bufferSize);
	private CharBuffer chBuf=CharBuffer.allocate(bufferSize);
	private Charset decoder= Charset.defaultCharset();
	private boolean hasNext;
	private String nextString,curString;
	private long append=0L,readBuffer=0L,decode=0L,search=0L,extractString=0L;
	private int start,stop;
	
	public LineReader(String file)throws Exception{
//		File file= new File(fileName);
		inChannel=(new FileInputStream(file)).getChannel();
//		System.err.println(readBuffer());
		readBuffer();
//		strBuf= new StringBuffer();
//		buf.flip();
		decode();
//		System.out.println(buf);
//		append();
//		buf.clear();
//		System.out.println(strBuf);
		nextString="";
		read();
	}
	
	public String readLine()throws Exception{
		if(!hasNext){
			throw new Exception("No more lines...");
		}
		curString=nextString;
		nextString="";
		read();
//		System.out.println(curString);
		return curString;
	}
	
	public boolean hasNext(){
		return hasNext;
	}
	
	private void append(){
		append-=System.currentTimeMillis();
//		while(buf.position()<buf.limit()){
//			strBuf.append((char)buf.get());
//		}
		//24747 and 5047
//		strBuf.append(decoder.decode(buf));
		//4377 and 4165
		chBuf.position(start);
		chBuf=chBuf.slice().append(decoder.decode(buf));
//		chBuf.flip();
		start=0;
//		chBuf.
		append+=System.currentTimeMillis();
	}
	
	private int readBuffer()throws Exception{
		readBuffer-=System.currentTimeMillis();
		buf.clear();
		final int ret=inChannel.read(buf);
		readBuffer+=System.currentTimeMillis();
		return ret;
	}
	
	private void decode(){
		decode-=System.currentTimeMillis();
		buf.flip();
		chBuf.clear();
		chBuf=decoder.decode(buf);
//		chBuf.flip();
		start=0;
		decode+=System.currentTimeMillis();
	}
	
	private int search(){
		search-=System.currentTimeMillis();
//		System.err.println(start);
//		System.err.println(chBuf.limit());
//		System.err.println(chBuf.toString());
		for(int i=0;i<chBuf.limit()-start;i++){
			final char c=chBuf.get();
//			System.err.println(""+c);
			if(c=='\n'){
				search+=System.currentTimeMillis();
				return i;
			}
		}
		search+=System.currentTimeMillis();
		return -1;
	}
	
	private void extractString(){
//		System.err.println(chBuf);
		extractString-=System.currentTimeMillis();
		chBuf.position(start);
		nextString+=chBuf.subSequence(0, stop).toString();
		//setup for next string
		start=start+stop+1;
//		start=start>chBuf.limit()?chBuf.limit():start;
//		chBuf.position(start);
		extractString+=System.currentTimeMillis();
	}
	
	public void printTimes(){
		System.err.println("LineReader:");
//		System.err.println("\tAppend: "+append);
		System.err.println("\tExtractString: "+extractString);
		System.err.println("\tSearch: "+search);
		System.err.println("\tDecode: "+decode);
		System.err.println("\tReadBuffer: "+readBuffer);
	}
	
	private void read()throws Exception{
		hasNext=false;
		stop=search();
//		System.err.println(stop);
		if(stop==-1){
			if(readBuffer()==-1){
				if(chBuf.limit()>0){
					stop=chBuf.limit()-start;
//					System.err.println(chBuf.limit());
					extractString();
					hasNext=true;
					chBuf.limit(0);
				}
			}else{
				stop=chBuf.limit()-start;
				extractString();
//				chBuf.position(start);
				decode();
				read();
			}
		}else{
//			System.err.println(chBuf.slice().toString());
			extractString();
			chBuf.position(start);
			hasNext=true;
		}
		
		
//		
////		boolean done=false;
//		int pos=strBuf.indexOf("\n");
////		System.err.println(pos);
//		if(pos==-1){
//			//read more
//			if(readBuffer()==-1){
//				if(buf.limit()>0){
//					System.err.println(buf.limit());
//					buf.flip();
////					strBuf.append(decoder.decode(buf));
////					strBuf.append(buf.asCharBuffer());
//					append();
////					System.out.println(strBuf.length());
//					buf.clear();
//					if(strBuf.length()>0){
////						read();
//					}
//				}else if(strBuf.length()>0){
//					//extract the last line
//					nextString=strBuf.toString();
//					hasNext=true;
//					strBuf.setLength(0);
//				}
//			}else{
//				buf.flip();
////				strBuf.append(decoder.decode(buf));
////				strBuf.append(buf.asCharBuffer());
//				append();
//				buf.clear();
//				read();
//			}
//		}else{
//			//extract substring
//			nextString=strBuf.substring(0, pos);
////			System.out.println(nextString);
//			strBuf.delete(0, pos+1);
//			hasNext=true;
//		}
		
//		while(!done){
//			if(buf.position()==buf.limit()){
////				System.err.println("Fill buffer..."+strBuf.length());
//				buf.clear();
//				if(inChannel.read(buf)==-1){
//					if(buf.limit()>0){
//						buf.flip();
////						read();
//						break;
//					}
//					if(strBuf.length()>0){
//						System.err.println("EOF");
//						nextString=strBuf.toString();
//						strBuf.setLength(0);
//						hasNext=true;
//					}else{
//						hasNext=false;
//					}
//					break;
//				}
//				buf.flip();
//			}
//			final char c=(char)buf.get();
//			if(c=='\n'){
//				nextString=strBuf.toString();
////				System.err.println(nextString);
//				strBuf.setLength(0);
//				hasNext=true;
//				done=true;
//			}else{
//				strBuf.append(c);
//			}
//		}
	}
}
