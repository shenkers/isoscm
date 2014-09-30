package executable;

import java.io.File;

import com.beust.jcommander.IStringConverter;

public class FileConverter implements IStringConverter<File> {

	public File convert(String value) {
		return new File(value).getAbsoluteFile();
	}

}
