name := "isoscm"

crossPaths := false

autoScalaLibrary := false

version := "2.0.5"

libraryDependencies += "org.xerial.snappy" % "snappy-java" % "1.1.1.7"

libraryDependencies += "com.github.samtools" % "htsjdk" % "1.130"

//libraryDependencies += "org.apache.commons" % "commons-lang3" % "3.4"

libraryDependencies += "commons-lang" % "commons-lang" % "2.6"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.4.1"

libraryDependencies += "commons-io" % "commons-io" % "2.4"

libraryDependencies += "com.beust" % "jcommander" % "1.29"

libraryDependencies += "org.apache.logging.log4j" % "log4j-core" % "2.2"

libraryDependencies += "org.apache.logging.log4j" % "log4j-api" % "2.2"

libraryDependencies += "com.beust" % "jcommander" % "1.29"

libraryDependencies += "org.apache.directory.studio" % "org.apache.commons.collections" % "3.2.1"

assemblyJarName in assembly := "IsoSCM.jar"

mainClass in assembly := Some("executable.IsoSCM")
