#ifndef DPGT_JNI_UTILS_HPP
#define DPGT_JNI_UTILS_HPP

#include <vector>
#include <string>
#include "jni.h"
#include "jni_md.h"


void JavaStringArrayToCppVec(JNIEnv *env, jobjectArray jarray, std::vector<std::string> &cppvec);


#endif  // DPGT_JNI_UTILS_HPP
