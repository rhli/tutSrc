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

JNIEXPORT jint JNICALL Java_JniTest_retValue (JNIEnv * a, jobject ab) {
  return 11;
}

JNIEXPORT void JNICALL Java_JniTest_printArr (JNIEnv * a, jobject ab, jintArray input) {
  jsize len = (*a) -> GetArrayLength(a, input);
  jint* arr = (*a) -> GetIntArrayElements(a, input, 0);
  jint i;
  printf("print array using JNI: \n");
  for (i = 0; i < len; i ++) {
    printf("%d ", arr[i]);
  }
  printf("\n");
}

JNIEXPORT jintArray JNICALL Java_JniTest_returnArr (JNIEnv * a, jobject ab, jint size) {
  jintArray arr;
  jintArray result = (*a) -> NewIntArray(a, size);
  jint i;
  if (!result) return NULL;
  jint fill[size];
  for (i = 0; i < size; i ++) {
    fill[i] = i + 1;
  }
  (*a) -> SetIntArrayRegion(a, result, 0, size, fill);
  return result;
}

JNIEXPORT void JNICALL Java_JniTest_square (JNIEnv * a, jobject ab, jobjectArray input, jobjectArray output) {
  jint inputLength = (*a) -> GetArrayLength(a, input);
  char** inputArr = (char**)calloc(sizeof(char**), inputLength); 
  char** outputArr = (char**)calloc(sizeof(char**), inputLength); 
  int i, j;
  for (i = 0; i < inputLength; i ++) {
    jobject jInput = (*a) -> GetObjectArrayElement(a, input, i);
    inputArr[i] = (*a) -> GetDirectBufferAddress(a, jInput);

    jobject jOutput = (*a) -> GetObjectArrayElement(a, output, i);
    outputArr[i] = (*a) -> GetDirectBufferAddress(a, jOutput);
  }

  for (i = 0; i < inputLength; i ++) {
    for (j = 0; j < inputLength; j ++) {
      outputArr[i][j] = inputArr[i][j] * inputArr[i][j];
    }
  }
}

