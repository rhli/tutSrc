package EARTest;
import EARTest.EARLayoutGen; 
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;

public class Main{
  public static void main(String[] args){
    Process p;
    //String command = new String("/bin/bash /home/rhli/playground/earTest/test.sh");
    String command = new String("touch abc");
    StringBuffer output=new StringBuffer();
    Map<String,String> envMap = System.getenv();
    System.out.println(envMap.get("HADOOP_HOME"));
    try {
      p = Runtime.getRuntime().exec(command);
      p.waitFor();
      //BufferedReader reader = 
      //  new BufferedReader(new InputStreamReader(p.getInputStream())); 
      //String line = ""; 
      //while ((line = reader.readLine())!= null) { 
      //  output.append(line + "\n"); 
      //}
      //System.out.println(output.toString());
    } catch(Exception e) {
    }
    //int repFac=2;
    //int ecK=10;
    //int ecN=14;
    //EARLayoutGen elg=new EARLayoutGen(ecK,ecN,repFac,12,1);
    //int[] output=new int[ecK*repFac];
    //System.out.println(System.nanoTime()/1000);
    //elg.SOPwoCoreRack(0,output);
    //System.out.println(System.nanoTime()/1000);
    //for(int i=0;i<10;i++){
    //  System.out.println("block " + i + ": " + output[i*repFac] + " " + output[i*repFac+1]);
    //}
  }
}
