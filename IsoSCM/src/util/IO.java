package util;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JWindow;

public class IO {
	
	public static class BroadcastPrintstream{
		Collection<PrintStream> streams;
		public BroadcastPrintstream(Collection<PrintStream> streams) {
			this.streams=streams;
		}
		
		public void print(String s){
			for(PrintStream p : streams){
				p.print(s);
			}
		}
		
		public void println(String s){
			for(PrintStream p : streams){
				p.print(s + "\n");
			}
		}
		
		public void println(){
			for(PrintStream p : streams){
				p.print("\n");
			}
		}

		public void close() {
			for(PrintStream p : streams){
				p.close();
			}
		}
	}

	//	public static void main(String[] args) throws IOException {
	//		DataOutputStream daos = new DataOutputStream(new FileOutputStream(new File("test")));
	//		daos.writeLong(1);
	//		daos.writeLong(3);
	//		daos.writeLong(7);
	//		daos.writeLong(8);
	//		daos.writeLong(100);
	//		daos.close();
	//		FileList fl = new FileList(new File("test"),0,5*(Long.SIZE/Byte.SIZE));
	//		System.out.println(binarySearch(fl, 100l));
	//	}

	public static int binarySearch(FileList sorted, long toFind) throws IOException{
		int length = sorted.size();

		int imin=0, imax=length-1;
		int imid = (imin + imax) / 2;

		// continue searching while [imin,imax] is not empty
		while (imax >= imin)
		{
			long val = sorted.getLong(imid);

			// determine which subarray to search
			if(val < toFind){
				// change min index to search upper subarray
				imin = imid + 1;
			}
			else if (val > toFind){
				// change max index to search lower subarray
				imax = imid - 1;
			}else{
				// key found at index imid
				return imid;
			}

			/* calculate the midpoint for roughly equal partition */
			imid = (imin + imax) / 2;
		}
		//		System.out.println(Util.list(imin,imid,imax));
		// key not found
		return -(imin + 1);
	}

	public static class FileList{
		File f;
		long start, end;
		int size;
		byte[] buffer;
		int longSize;
		int byteSize;

		public FileList(File f, long start, long end, int byteSize) throws FileNotFoundException {
			this.f = f;
			this.start = start;
			this.end = end;
			this.byteSize = byteSize;
			longSize = Long.SIZE/Byte.SIZE;
			size = (int) ((end-start)/(byteSize));
			buffer = new byte[byteSize];
		}

		public Long getLong(int index) throws IOException {
			//			System.out.println(index);
			RandomAccessFile raf = new RandomAccessFile(f, "r");
			raf.seek(index*byteSize + start);
			long l = raf.readLong();
			raf.close();
			return l;
		}

		public byte[] get(int index) throws IOException {
			//			System.out.println(index);
			RandomAccessFile raf = new RandomAccessFile(f, "r");
			raf.seek(index*byteSize + start);
				raf.seek(index*byteSize + start);
			raf.readFully(buffer);
			raf.close();
			return buffer;
		}

		public int size() {
			return size;
		}
	}

	public static void saveImage(Component component, String fileName, String formatName, int width, int height) throws IOException{
		JWindow jw = new JWindow();
		jw.add(component);
		jw.pack();
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_4BYTE_ABGR);
		Graphics g = img.getGraphics();
		component.paint(img.getGraphics());
		g.dispose();
		jw.dispose();
		ImageIO.write(img, formatName, new File(fileName));
	}

	public static class LineTokenizer implements Iterable<String[]>, Iterator<String[]>{

		Scanner scanner;
		String regex;

		public LineTokenizer(Scanner scanner, String regex) {
			this.scanner = scanner;
			this.regex = regex;
		}

		public Iterator<String[]> iterator() {
			return this;
		}

		public boolean hasNext() {
			return scanner.hasNextLine();
		}

		public String[] next() {
			return scanner.nextLine().split(regex,-1);
		}

		public void remove() {

		}

	}

	public static Scanner bufferedScanner(String fileName) throws FileNotFoundException{
		return new Scanner(new BufferedInputStream(new FileInputStream(fileName)));
	}

	public static Scanner bufferedScanner(File file) throws FileNotFoundException{
		return new Scanner(new BufferedInputStream(new FileInputStream(file)));
	}
	
	public static Scanner bufferedScanner(InputStream source) {
		return new Scanner(new BufferedInputStream(source));
	}
	
	public static FileOutputStream OutFileOutputStream(String fileName, boolean overwrite) throws FileNotFoundException{
		File f = new File(fileName);
		if(!overwrite && f.exists())
			throw new IllegalArgumentException("File exists");
		return new FileOutputStream(fileName);
	}
	
	public static PrintStream OutPrintstream(String fileName, boolean overwrite) throws FileNotFoundException{
		File f = new File(fileName);
		if(!overwrite && f.exists())
			throw new IllegalArgumentException("File exists");
		return new PrintStream(f);
	}
	
	public static PrintStream OutPrintstream(File f, boolean overwrite) throws FileNotFoundException{
		if(!overwrite && f.exists())
			throw new IllegalArgumentException("File exists");
		return new PrintStream(f);
	}

	public static PrintStream OutPrintstream(File infile, String suffix, boolean overwrite) throws FileNotFoundException{
		File f = new File(infile.getAbsolutePath()+suffix);
		if(!overwrite && f.exists())
			throw new IllegalArgumentException("File exists");
		return new PrintStream(f);
	}

	public static PrintStream OutPrintstream(String infileName, String suffix, boolean overwrite) throws FileNotFoundException{
		File f = new File(infileName+suffix);
		if(!overwrite && f.exists())
			throw new IllegalArgumentException(Util.sprintf("File \"%s\" exists", f.getAbsolutePath()));
		return new PrintStream(f);
	}

	public static PrintStream bufferedPrintstream(String fileName) throws FileNotFoundException{
		return new PrintStream(new BufferedOutputStream(new FileOutputStream(fileName)));
	}

	public static PrintStream bufferedPrintstream(File file) throws FileNotFoundException{
		return new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));
	}

	public static class GetOpts{
		Map<String,String> opts;
		public GetOpts(String[] args){
			opts = new HashMap<String,String>();
			String key = null;
			Pattern p = Pattern.compile("^-");

			for(int i=0; i<args.length; i++){
				Matcher m = p.matcher(args[i]); // get a matcher object
				if(m.find()){
					key = args[i].substring(1);	
					opts.put(key,null);
				}
				else{
					opts.put(key,args[i]);
				}
			}
		}

		public String get(String option){
			if(!opts.containsKey(option))
				throw new IllegalArgumentException("Command line option '-" + option + "' is required");
			return opts.get(option);
		}

		public <T> T get(String option, T type) throws IllegalArgumentException, InstantiationException, IllegalAccessException, InvocationTargetException, SecurityException, NoSuchMethodException{

			if(!opts.containsKey(option))
				throw new IllegalArgumentException("Command line option '-" + option + "' is required");
			Constructor<T> c = (Constructor<T>) type.getClass().getConstructor(String.class);

			return c.newInstance(opts.get(option));
		}

		public boolean contains(String option){
			return opts.containsKey(option);
		}

		public <T> T optional(String option, T alternative) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException {
			if(!opts.containsKey(option))
				return alternative;
			else
				return get(option, alternative);
		}
	}
	
	public static JFrame frame(Component c){
		JFrame f = new JFrame();
		f.getContentPane().add(c);
		return f;
	}
	
	public static void display(JFrame c){
		c.pack();
		c.setResizable(true);
		c.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		c.setLocationRelativeTo(null);
		c.setVisible(true);
	}


	public static JFrame show(Component c){
		JFrame f = new JFrame();
		f.getContentPane().add(c);
		f.pack();
		f.setResizable(true);
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		f.setLocationRelativeTo(null);
		f.setVisible(true);
		return f;
	}
	
	public static JFrame show(Component c, boolean resizable){
		JFrame f = new JFrame();
		f.getContentPane().add(c);
		f.pack();
		f.setResizable(resizable);
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		f.setLocationRelativeTo(null);
		f.setVisible(true);
		return f;
	}

	public static void shutdownButton(){
		JFrame f = new JFrame();
		JButton b = new JButton("Terminate");
		b.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				Runtime.getRuntime().exit(0);
			}
		});
		f.getContentPane().add(b);
		f.pack();
		f.setResizable(false);
		f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		f.setLocationRelativeTo(null);
		f.setVisible(true);
	}

	public static abstract class ComponentAdaptor extends Component{

		private static final long serialVersionUID = 1L;

		public abstract void paint(Graphics g);

	}

	public static void setPrefferedSize(Component c, int w, int h) {
		c.setPreferredSize(new Dimension(w,h));
	}

	public static Dimension dimension(int width, int height) {
		return new Dimension(width, height);
	}

	
}


