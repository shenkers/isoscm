package executable;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;


public class ConfigurationIO {
	
	/**
	 * 
	 * @param command
	 * @param configuration_xml .xml file to which the configuration will be written
	 * @throws JAXBException
	 */
	public static void writeConfiguration(AssembleCommand command, File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(AssembleCommand.class);
		
		Marshaller m = jaxbc.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

		m.marshal(command, configuration_xml);
	}
	
	/**
	 * reads the command parameters from the configuration .xml file
	 * @param command
	 * @param configuration_xml
	 * @throws JAXBException
	 */
	public static AssembleCommand readAssemblyConfiguration( File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(AssembleCommand.class);
		
		Unmarshaller um = jaxbc.createUnmarshaller();
		
		return (AssembleCommand) um.unmarshal(configuration_xml);
	}
	
	/**
	 * 
	 * @param command
	 * @param configuration_xml .xml file to which the configuration will be written
	 * @throws JAXBException
	 */
	public static void writeConfiguration(CompareCommand command, File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(CompareCommand.class);
		
		Marshaller m = jaxbc.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

		m.marshal(command, configuration_xml);
	}
	
	/**
	 * reads the command parameters from the configuration .xml file
	 * @param command
	 * @param configuration_xml
	 * @throws JAXBException
	 */
	public static CompareCommand readCompareConfiguration(File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(CompareCommand.class);
		
		Unmarshaller um = jaxbc.createUnmarshaller();
		
		return (CompareCommand) um.unmarshal(configuration_xml);
	}
	
	/**
	 * 
	 * @param command
	 * @param configuration_xml .xml file to which the configuration will be written
	 * @throws JAXBException
	 */
	public static void writeConfiguration(EnumerateCommand command, File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(EnumerateCommand.class);
		
		Marshaller m = jaxbc.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

		m.marshal(command, configuration_xml);
	}
	
	/**
	 * reads the command parameters from the configuration .xml file
	 * @param command
	 * @param configuration_xml
	 * @throws JAXBException
	 */
	public static EnumerateCommand readEnumerateConfiguration(File configuration_xml) throws JAXBException{
		JAXBContext jaxbc = JAXBContext.newInstance(AssembleCommand.class);
		
		Unmarshaller um = jaxbc.createUnmarshaller();
		
		return (EnumerateCommand) um.unmarshal(configuration_xml);
	}
	
	@XmlRootElement
	static class A{
		int age;
		String name;
		double weight;
		File f;
		
		public A() {
		
		}

		public A(int age, String name, double weight, File f) {
			this.age=age;
			this.name=name;
			this.weight=weight;
			this.f = f;
		}
		
		public String getName(){
			return name;
		}
		
		public double getWeight(){
			return weight;
		}
		
		@XmlElement
		public void setName(String name){
			this.name=name;
		}
		
		@XmlElement
		public void setWeight(double weight){
			this.weight=weight;
		}
		
		public File getFile(){
			return f;
		}
		
		@XmlElement
		public void setFile(File f){
			this.f=f;
		}
	}
	
	@XmlRootElement(name="list")
	static class cA{
		List<A> a;
		
		public cA() {
			a=new LinkedList<A>();
		}
		
		@XmlAnyElement(lax=true)
		public List<A> getA(){
			return a;
		}
	}

	public static void main(String[] args) throws JAXBException, IOException {
		
		JAXBContext jc = JAXBContext.newInstance(cA.class,A.class);
		Marshaller m = jc.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

		cA c = new cA();
		c.getA().add(new A(13, "name a", 10.2,new File("a")));
		c.getA().add(new A(13, "name b", 10.2,new File("b")));
//		File f = new File("test.xml");
//		m.marshal(c, f);
		m.marshal(c, System.out);
		Unmarshaller um = jc.createUnmarshaller();
//		cA cu = (cA) um.unmarshal(f);
//		for(A a : cu.getA()){
//			System.out.printf("unmarshed %s\n", a.name);
//		}
	}
}
