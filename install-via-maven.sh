#/usr/bin/bash

root_dir=$(readlink -f "$(dirname $0)")
echo $root_dir 

#Cleaning
echo "Cleaning..."
mvn clean

#Adding phaseIIFilter
echo -e "\n\nAdding Phase II Filter JAR..."
mvn org.apache.maven.plugins:maven-install-plugin:3.0.0-M1:install-file -Dfile="${root_dir}/lib/phase2filter-1.0.2.jar" -DgroupId=djoumbou -DartifactId=phase2filter -Dversion=1.0.2 -Dpackaging=jar

#Adding Cypreact
echo -e "\n\nAdding CypReact JAR..."
mvn org.apache.maven.plugins:maven-install-plugin:3.0.0-M1:install-file -Dfile="${root_dir}/lib/CypReact.jar" -DgroupId=wishartlab -DartifactId=CypReact -Dversion=1.0.0-SNAPSHOT -Dpackaging=jar

#Adding CyProduct
echo -e "\n\nAdding CyProduct JAR..."
mvn org.apache.maven.plugins:maven-install-plugin:3.0.0-M1:install-file -Dfile="${root_dir}/lib/CyProduct.jar" -DgroupId=wishartlab -DartifactId=CyProduct -Dversion=1.0.0-SNAPSHOT -Dpackaging=jar

#Building package
echo -e "\n\nBuilding Package..." 
mvn package

#Moving jar file
echo -e "Moving .jar file to the root directory ..."
mv target/*.jar ./

echo "Done!"


