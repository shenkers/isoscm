package executable;

import java.io.File;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.beust.jcommander.Parameter;

@XmlRootElement(name="CompareCommand")
@XmlAccessorType(XmlAccessType.PROPERTY)
public class CompareCommand{	


	@Parameter(names="-x1", description="XML configuration file for the assembly step from sample 1", converter=FileConverter.class)
	File assemblyXml1;

	@Parameter(names="-x2", description="XML configuration file for the assembly step from sample 2", converter=FileConverter.class)
	File assemblyXml2;

	@Parameter(names="-base", description="output files will be written to [dir]/[base].txt and [dir]/[base].gtf")
	String base;

	@Parameter(names="-dir", description="The output directory", converter=FileConverter.class)
	public File dir;

	public File getAssemblyXml1() {
		return assemblyXml1;
	}

	public File getAssemblyXml2() {
		return assemblyXml2;
	}
	
	public String getBase() {
		return base;
	}	

	public File getDir() {
		return dir;
	}

	@XmlElement
	public void setAssemblyXml1(File assemblyXml1) {
		this.assemblyXml1 = assemblyXml1;
	}

	@XmlElement
	public void setAssemblyXml2(File assemblyXml2) {
		this.assemblyXml2 = assemblyXml2;
	}

	@XmlElement
	public void setBase(String base) {
		this.base = base;
	}

	@XmlElement
	public void setDir(File dir) {
		this.dir = dir;
	}
}

