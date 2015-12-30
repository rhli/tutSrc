/**
 * @file NativeStair.c
 * @brief Native C implementation of stair codes.
 * @author Mingqiang Li, mingqiangli.cn@gmail.com
 * @version 0.1
 * @date 2014-12-04
 */
#include "NativeStair.h"

int n_original;
int r_original;
int m_original;
int m_partial_original;
int *error_vector_original = NULL;

bool *parity_layout_original = NULL;

int word_size_in_gf;
gf_t gf_obj;

int *row_coding_matrix_original = NULL;
int *column_coding_matrix_original = NULL;
/* Added by RH begins */
int total_parity_count;
int* stripeDataMap;
int* stripeParityMap;
int* parityStripeMap;
int* dataStripeMap;
/* Added by RH ends */

int best_encoding_method;

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeInit - initialize all related variables
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param stripeSize - number of data symbols in a stripe
 * @param paritySize - number of parity symbols in a stripe
 * @param numLocalGroup - number of local groups
 * @param sizeLocalGroup - size of each local group
 * @param localParitySize - number of local parity symbols in each local group
 * @param globalParityVector - global parity vector
 * @param globalParityVectorLen - length of global parity vector
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool nativeInit(int stripeSize, int paritySize, int numLocalGroup, int sizeLocalGroup, 
 int localParitySize, int* globalParityVector, int globalParityVectorLen) {
  fprintf(stderr,"nativeInit(): checkpoint 1\n");
	int *global_parity_vector = globalParityVector;
	int global_parity_vector_size = globalParityVectorLen;
	int i, j;

	n_original = sizeLocalGroup;
	r_original = numLocalGroup;
	m_original = localParitySize;
  /* Added by RH Dec 20th begins */
  total_parity_count = m_original*r_original;
  /* Added by RH Dec 20th ends */
	m_partial_original = global_parity_vector[global_parity_vector_size - 1];

	/*check n_original, r_original, m_original, and m_partial_original configurations*/
	if (n_original <= 0) {
		fprintf(stderr, "Error: bad n_original (%d), which should be > 0!\n", n_original);
		return false;
	}
	if (r_original <= 0) {
		fprintf(stderr, "Error: bad r_original (%d), which should be > 0!\n", r_original);
		return false;
	}
	if ((m_original <= 0) || (m_original >= n_original)) {
		fprintf(stderr, "Error: bad m_original (%d), which should be in (0, %d)!\n", m_original, n_original);
		return false;
	}	
	if ((m_partial_original <= 0) || (m_partial_original > n_original - m_original)) {
		fprintf(stderr, "Error: bad m_partial_original (%d), which should be in (0, %d]!\n", 
				m_partial_original, n_original - m_original);
		return false;
	}	
  //fprintf(stderr,"nativeInit(): checkpoint 2\n");

	/*check global_parity_vector*/
	int pre = 0;
	for (i = 0; i < global_parity_vector_size; i++) {	
		if ((global_parity_vector[i] <= 0) || (global_parity_vector[i] >= n_original - m_original)) {
			fprintf(stderr, "Error: bad %d-th element (%d) in the input extra vector, which should be in (0, %d)!\n", 
					i, global_parity_vector[i], n_original - m_original);
			return false;
		}
		if (global_parity_vector[i] < pre) {
			fprintf(stderr, "Error: the %d-th element (%d) in the input extra vector is not in \
					monotonically increasing order!\n", i, global_parity_vector[i]);
			return false;
		}

		pre = global_parity_vector[i];
	}

	/*generate error_vector_original*/
	error_vector_original = TypeAlloc(int, m_partial_original);
	for (i = 0; i < m_partial_original; i++) {	
		error_vector_original[i] = 0;
		j = global_parity_vector_size - 1;

		while ((j >= 0) && (global_parity_vector[j] >= m_partial_original - i)) {
			error_vector_original[i]++;
			j--;			
		}
    total_parity_count+=error_vector_original[i];
    printf("error_vector_original[%d]=%d\n",i,error_vector_original[i]);
	}

	InitParityLayout(n_original, r_original, m_original, m_partial_original,
			error_vector_original, &parity_layout_original,
      &parityStripeMap,&stripeParityMap,&dataStripeMap,&stripeDataMap);

	if(!InitGF(n_original, r_original, m_original, m_partial_original, 
				error_vector_original, &word_size_in_gf, &gf_obj)) {
		return false;
	}

	InitCodingMatrices(gf_obj, n_original, r_original, m_original, m_partial_original,
			&error_vector_original, &row_coding_matrix_original, &column_coding_matrix_original);
  if(row_coding_matrix_original==NULL){
    fprintf(stderr,"NativeStair.c: row_coding_matrix_original not initialized!!\n");
  }

	best_encoding_method = ChooseEncodingMethod(n_original, r_original,
			m_original, m_partial_original, error_vector_original, parity_layout_original) ;

	fprintf(stderr, "\nSucceed in initializing a transposed stair code with the following parameters:-)\n");
	fprintf(stderr, "Parameters:\n");
	fprintf(stderr, "      n_original: %d\n", n_original);			
	fprintf(stderr, "      r_original: %d\n", r_original);							
	fprintf(stderr, "      m_original: %d\n", m_original);					
	fprintf(stderr, "      error_vector_original: (");
	for (i = 0; i < m_partial_original - 1; i++) {		
		fprintf(stderr, "%d, ", error_vector_original[i]);	
	}					
	fprintf(stderr, "%d)\n", error_vector_original[m_partial_original - 1]);	
	fprintf(stderr, "      word_size_in_gf: %d\n", word_size_in_gf);					
	fprintf(stderr, "      row_coding_matrix_original: (see below)\n");				
	for (i = 0; i < m_original + m_partial_original; i++) {			
		fprintf(stderr, "         | ");			
		for (j = 0; j < n_original - m_original; j++) {				
			fprintf(stderr, "%3d ", row_coding_matrix_original[(n_original - m_original) * i + j]);			
		}			
		fprintf(stderr, "|\n");		
	}							
	fprintf(stderr, "      column_coding_matrix_original: (see below)\n");				
	for (i = 0; i < error_vector_original[m_partial_original - 1]; i++) {			
		fprintf(stderr, "         | ");			
		for (j = 0; j < r_original; j++) {				
			fprintf(stderr, "%3d ", column_coding_matrix_original[r_original * i + j]);			
		}			
		fprintf(stderr, "|\n");
	}	
	if (best_encoding_method == STD_METHOD) {
		fprintf(stderr, "      best_encoding_method: STD_METHOD\n");	
	}

	if (best_encoding_method == UP_METHOD) {
		fprintf(stderr, "      best_encoding_method: UP_METHOD\n");	
	}

	if (best_encoding_method == DOWN_METHOD) {
		fprintf(stderr, "      best_encoding_method: DOWN_METHOD\n");	
	}			
	fprintf(stderr, "      decoding method: DOWN_METHOD\n");		
	fprintf(stderr, "\n");	

	return true;
}

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeCleanup - clean up all allocated space
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 */
//JNIEXPORT void JNICALL Java_org_apache_hadoop_util_NativeStair_nativeCleanup
//(JNIEnv *env, jclass clazz) {
//	free(error_vector_original);
//	free(parity_layout_original);
//
//	/*free the gf_t object*/		
//	gf_free(&gf_obj, 1);
//
//	free(row_coding_matrix_original);
//	free(column_coding_matrix_original);
//
//	fprintf(stderr, "\nWe have cleaned up all space allocated for the stair code!\n");
//	fprintf(stderr, "\n");	
//}
//
///**
// * @brief Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize - get symbol size
// *
// * @param env - JNI interface pointer
// * @param clazz - Java class object
// *
// * @return - symbol size
// */
//JNIEXPORT int JNICALL Java_org_apache_hadoop_util_NativeStair_nativeSymbolSize
//(JNIEnv *env, jclass clazz) {
//	return word_size_in_gf;
//}
//
/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeEncodeBulk - encode data into parity in bulk
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param inputBuffers - input buffers that store data 
 * @param numInputBuffers - number of input buffers
 * @param outputBuffers - output buffers that store parity 
 * @param numOutputBuffers - number of output buffers
 * @param dataLength - length of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool nativeEncodeBulk (char** data, int numInputBuffers, 
 char** coding, int numOutputBuffers, int dataLength) {
	int i;

	if ((dataLength <= 0) || (dataLength % (word_size_in_gf/8) != 0)) {
		fprintf(stderr, "Error: bad dataLength (%d), which should be a positive multiple of %d!\n", 
				dataLength, word_size_in_gf/8);
		return false;
	}
	//return true;

	if(!EncodeBulk(gf_obj, n_original, r_original, m_original, m_partial_original, 
				error_vector_original, parity_layout_original, best_encoding_method, 
				row_coding_matrix_original, column_coding_matrix_original, 
				data, coding, dataLength)) {
		return false;
	}

	return true;
}

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode - determine locations 
 *                 to be read for decoding
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param erasedLocationsBitmap - bitmap of erased locations
 * @param locationsToReadBitmap - bitmap of locations to be read
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool nativeLocationsToReadForDecode (bool* erased_layout, int* read_layout,int locBitmapSize) {
//JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_util_NativeStair_nativeLocationsToReadForDecode
//(JNIEnv *env, jclass clazz, jintArray erasedLocationsBitmap, jintArray locationsToReadBitmap) {
  int i;
  bool* stripeErasedLayout=(bool*)calloc(sizeof(bool),n_original*r_original);
  for(i=0;i<n_original*r_original;i++){
    if(erased_layout[i]){
      if(i<total_parity_count){
        stripeErasedLayout[parityStripeMap[i]]=true;
      }else{
        stripeErasedLayout[dataStripeMap[i-total_parity_count]]=true;
      }
    }
  }
  for(i=0;i<n_original*r_original;i++){
    printf(stripeErasedLayout[i]?"stripeErasedLayout[%d]=1\n":"stripeErasedLayout[%d]=0\n",i);
  }
  int* stripeReadLayout=(int*)calloc(sizeof(int),n_original*r_original);
	if (!FindLocationsToReadForDecode(n_original, r_original, m_original, 
				m_partial_original, error_vector_original, stripeErasedLayout, stripeReadLayout)) {
    return false;
	}
  for(i=0;i<n_original*r_original;i++){
    if(stripeReadLayout[i]==1){
      printf("stripeReadLayout[%d]=1\n",i);
      if(stripeDataMap[i]!=-1){
        printf("stripeDataMap[%d]=%d\n",i,stripeDataMap[i]);
        read_layout[stripeDataMap[i]+total_parity_count]=1;
      }else{
        printf("stripeParityMap[%d]=%d\n",i,stripeParityMap[i]);
        read_layout[stripeParityMap[i]]=1;
      }
    }
  }
  return true;
}

/**
 * @brief Java_org_apache_hadoop_util_NativeStair_nativeDecodeBulk - decode erased blocks in bulk
 *
 * @param env - JNI interface pointer
 * @param clazz - Java class object
 * @param inputBuffers - input buffers that store the whole stripe 
 * @param numInputBuffers - number of input buffers
 * @param outputBuffers - output buffers that store the decoded blocks 
 * @param numOutputBuffers - number of output buffers
 * @param erasedLocationsArray - an array storing erased locations
 * @param numErasedLocations - number of erased locations
 * @param dataLength - length of each data/parity block
 *
 * @return - a boolean value to indicate if the operation is successful
 */
bool nativeDecodeBulk (char** stripe, int numInputBuffers, char** outputBufs,
 int numOutputBuffers, bool* erased_layout, int numErasedLocations, int dataLength){
	int i;

	if ((dataLength <= 0) || (dataLength % (word_size_in_gf/8) != 0)) {
		fprintf(stderr, "Error: bad dataLength (%d), which should be a positive multiple of %d!\n", 
				dataLength, word_size_in_gf/8);
		return false;
	}

  char** stripeLayout=(char**)calloc(sizeof(char*),n_original*r_original);
  for(i=0;i<total_parity_count;i++){
    stripeLayout[parityStripeMap[i]]=stripe[i];
  }
  for(i=total_parity_count;i<n_original*r_original;i++){
    stripeLayout[dataStripeMap[i-total_parity_count]]=stripe[i];
  }
  bool* stripeErasedLayout=(bool*)calloc(sizeof(bool),n_original*r_original);
  for(i=0;i<n_original*r_original;i++){
    if(erased_layout[i]){
      if(i<total_parity_count){
        stripeErasedLayout[parityStripeMap[i]]=true;
      }else{
        stripeErasedLayout[dataStripeMap[i-total_parity_count]]=true;
      }
    }
  }
  //fprintf(stderr,"before DecodeBulk()\n");
	if (!DecodeBulk(gf_obj, n_original, r_original, m_original, m_partial_original, 
				error_vector_original, parity_layout_original, row_coding_matrix_original, 
				column_coding_matrix_original, stripeLayout, dataLength, stripeErasedLayout, outputBufs)) {
		return false;
	}

	return true;
}
