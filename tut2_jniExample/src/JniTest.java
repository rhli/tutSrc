import java.nio.ByteBuffer;

public class JniTest{ 
  static{ 
    System.loadLibrary("jnitest"); 
  }

  private native int retValue();
  private native void sayHello();
  private native void passArg(int arg);
  private static native void printArr(int[] arr);
  private static native int[] returnArr(int size);
  private static native void square(ByteBuffer[] input, ByteBuffer[] output);

  public static void main(String[] args){
    JniTest jt = new JniTest();
    // return primitive values
    System.out.println("retValue() returns: " + jt.retValue());

    // return array
    int[] retArr = returnArr(17);
    for (int i = 0; i < retArr.length; i ++) {
      System.out.println("retArr[" + i + "]=" + retArr[i]);
    }

    // hello world function
    jt.sayHello();
    
    // passing in argument
    jt.passArg(10);

    // passing in array
    int arr[] = {3,2,1};
    printArr(arr);

    // two dimensional array and direct access
    byte[][] input = {{1, 2}, {3, 4}};
    byte[][] output = {{0, 0}, {0, 0}};
    ByteBuffer[] inputBufs = new ByteBuffer[input.length];
    ByteBuffer[] outputBufs = new ByteBuffer[output.length];
    for (int i = 0; i < input.length; i ++) {
      inputBufs[i] = ByteBuffer.allocateDirect(input[i].length);
      inputBufs[i].put(input[i], 0, input[i].length);
    }
    for (int i = 0; i < output.length; i ++) {
      outputBufs[i] = ByteBuffer.allocateDirect(output[i].length);
    }
    square(inputBufs, outputBufs);
    for (int i = 0; i < output.length; i ++) {
      outputBufs[i].get(output[i]);
    }
    System.out.println("Square results:");
    for (int i = 0; i < output.length; i ++) {
      System.out.println(output[i][0] + " " + output[i][1]);
    }
  }
}
