import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
 
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.BlockLocation;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hdfs.DistributedFileSystem;
import org.apache.hadoop.hdfs.protocol.DatanodeInfo;

public class example{
  private static Configuration conf = new Configuration();
  private static FileSystem fileSystem;
  
  public static void mkdir(String dir) throws IOException{
    Path path = new Path(dir);
    if (fileSystem.exists(path)) {
      System.out.println("Dir " + dir + " already exists!");
      return;
    }
    fileSystem.mkdirs(path);
  }

  public static void main(String[] args) throws IOException {
    System.out.println("test");
    conf.addResource(new Path("/home/Projects/previous/hadoop-20-master/conf/core-site.xml"));
    conf.addResource(new Path("/home/Projects/previous/hadoop-20-master/conf/hdfs-site.xml"));
    conf.addResource(new Path("/home/Projects/previous/hadoop-20-master/conf/mapred-site.xml"));
    String dir = new String("apiTest");

    fileSystem = FileSystem.get(conf);
    mkdir(dir);
  }
}
