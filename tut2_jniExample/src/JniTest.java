public class JniTest{ 
  static{ 
    System.loadLibrary("jnitest"); 
  }

  private native void sayHello();
  private native void passArg(int arg);

  public static void main(String[] args){
    JniTest jt = new JniTest();
    jt.sayHello();
    jt.passArg(10);
  }
}
