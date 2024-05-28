#include "jni_utils.hpp"

void JavaStringArrayToCppVec(JNIEnv *env, jobjectArray jarray, std::vector<std::string> &cppvec)
{
    int n = env->GetArrayLength(jarray);
    cppvec = std::vector<std::string>(n, "");
    for (int i = 0; i < n; ++i) {
        jstring s1 = (jstring) (env->GetObjectArrayElement(jarray, i));
        const char *s2 = env->GetStringUTFChars(s1, NULL);
        cppvec[i] = s2;
        env->ReleaseStringUTFChars(s1, s2);
    }
}