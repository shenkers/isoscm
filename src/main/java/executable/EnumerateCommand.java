package executable;

import java.io.File;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.IntegerConverter;

@XmlRootElement(name="EnumerateCommand")
@XmlAccessorType(XmlAccessType.PROPERTY)
public class EnumerateCommand{	
    
	@Parameter(names="-x", description="configuration XML file from the assembly step", converter=FileConverter.class, required=true)
	File assemblyXml;
	
	@Parameter(names="-max_isoforms", description="loci with more than this number of isoforms will be skipped", converter=IntegerConverter.class, required=true)
	Integer max_paths;
	
	@XmlElement
	public void setAssemblyXml(File assemblyXml) {
		this.assemblyXml = assemblyXml;
	}

	@XmlElement
	public void setMax_paths(Integer max_paths) {
		this.max_paths = max_paths;
	}

	public File getAssemblyXml() {
		return assemblyXml;
	}
	
	public Integer getMax_paths() {
		return max_paths;
	}
}
