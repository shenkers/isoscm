package executable;

import java.io.File;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.beust.jcommander.Parameter;

@XmlRootElement(name="CompareCommand")
@XmlAccessorType(XmlAccessType.PROPERTY)
public class DiffCommand{	


	@Parameter(names="-x", description="XML configuration file from the 'compare' step", converter=FileConverter.class, required=true)
	File compareXml;

	@Parameter(names="-G", description="Ensembl format GTF to which the assembled exons will be compared", required=true)
	String refGtf;

	public File getAssemblyXml() {
		return compareXml;
	}

	@XmlElement
	public void setAssemblyXml(File assemblyXml) {
		this.compareXml = assemblyXml;
	}


}

