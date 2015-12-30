#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include "JniTest.h"


JNIEXPORT void JNICALL Java_JniTest_sayHello (JNIEnv * a, jobject ab) { 
  printf("This is a JNI program\n");
}

JNIEXPORT void JNICALL Java_JniTest_passArg (JNIEnv * a, jobject ab, jint arg) { 
  printf("received arg is %d\n", arg);
}

