//
// rjava.c:  rlab's interface to Java Virtual Machine (JVM)
// Marijan Kostrun, VIII-2012
//
// This file is a part of RLaB + rlabplus
// Copyright (C) 2012  Marijan Kostrun
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ../COPYING

#include "rlab.h"
#include "listnode.h"
#include "ent.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "class.h"
#include "rfileio.h"
#include "mathl.h"

#ifdef HAVE_JVM  /* Set in config.h */

#include <jni.h>
static JavaVM *jvm = 0;
static JNIEnv *env;

#ifdef _WIN32
# define PATH_SEPARATOR ';'
#else /* UNIX */
# define PATH_SEPARATOR ':'
#endif

#define USER_CLASSPATH "." /* where Prog.class is */

#define MAXLEN_JAVA_DESC 2048
static char method_desc_rlab[MAXLEN_JAVA_DESC] = {'\0'};

/*******************************************************************************
 * Create a Java Virtual Machine as part of the RLaB Process...
 *
 * Args:  arguments to the Java VM
 * Examples:
 * j_opts = ["-Djava.compiler=NONE", "-Djava.class.path=../../myclasses", ...
 *  "-Djava.library.path=../../mylibs",  "-verbose:jni" ];
 * jvm.create( j_opts )
 ******************************************************************************/

Ent *
ent_jvm_CreateJVM (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDS *ropt=0;

  JavaVMOption *options=0;
  JavaVMInitArgs vm_args;

  int i, nopts=0, retval=1;

  // if JVM exists just return 1
  if (jvm)
  {
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar ((double) retval);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  // The user can put JVM options in the argument list to
  // createJVM ()
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
    {
      ropt  = ent_data(e1);
      nopts = ropt->nrow * ropt->ncol;
    }
  }

  // did user provide any options?
  if (nopts)
  {
    options = (JavaVMOption *) GC_MALLOC (sizeof (JavaVMOption) * nopts);
    for (i=0; i<nopts; i++)
      options[i].optionString = MdsV0(ropt, i);
  }

  // initialize JVM
  vm_args.version  = JNI_VERSION_1_2;
  vm_args.options  = options;
  vm_args.nOptions = nopts;
  vm_args.ignoreUnrecognized = JNI_TRUE;
  retval = JNI_CreateJavaVM(&jvm, (void **)&env, &vm_args);

  // start cleaning up
  ent_Clean(e1);
  if (options)
    GC_FREE (options);

  if (retval < 0)
    rerror ("ent_jvm_CreateJVM: failed to create Java VM !");

  return ent_Create_Rlab_Double((double) retval);
}

/*******************************************************************************
 * Destroy a Java Virtual Machine as part of the RLaB Process...
 *
 * Args:  does not take any.
 * Note: though this is supposed to destroy it, once destroyed, another
 *       cannot be started, or so it seems.
 ******************************************************************************/

Ent *
ent_jvm_DestroyJVM (int nargs, Datum args[])
{
  int retval=0;

  if (jvm)
  {
    (*jvm)->DestroyJavaVM(jvm);
    retval = 1;

    printf("Note: one cannot unload the JVM in JDK release 1.1 or Java 2 SDK release 1.2.\n");
    printf("DestroyJVM always returns an error code in these releses.\n");
    printf("(cf. The Java Native Interface Programmer's guide and Specification, p.86)\n");
  }

  return ent_Create_Rlab_Double((double) retval);
}

/*******************************************************************************
 * Start a Java appp. soon this will be multiple real functions. For now
 * its hardcoded to start something specific up.
 ******************************************************************************/

Ent *
ent_jvm_CallMethod (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0, *es=0;
  char *meth0=0, *method=0, *method_name=0, *method_class=0;
  char *method_desc=0, *desc_rval=0;
  MDS *s=0, *sval=0;
  MDR *r=0, *rval=0;
  int i, j, iretval=0, hn=0, nr, ns, ist_meth=0;
  double retval=0;
  ListNode *node=0;

  jclass jCls;
  jmethodID jMid;
  jobjectArray jObj=0, jRetarray=0;
  jvalue *jArgs=0;

  method_desc_rlab[0] = '\0';

  if (!nargs)
    rerror ("ent_jvm_CallMethod: at least one argument required");

  //
  // first argument is either:
  // 2. <<method;descriptor>> where
  //  'method' is the name of the method
  //  while 'descriptor' is JNI-formatted string containing types and order of the input
  //  arguments and the return type of the method
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_STRING)
    meth0 = class_char_pointer (e1);
  else if (ent_type (e1) == BTREE)
  {
    // arg1.method
    node = btree_FindNode (ent_data (e1), "method");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_STRING)
        meth0 = class_char_pointer(var_ent (node));
    }
    // arg1.descriptor
    node = btree_FindNode (ent_data (e1), "descriptor");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_STRING)
        method_desc = class_char_pointer(var_ent (node));
    }
  }

  if (!meth0)
    rerror ("ent_jvm_CallMethod: first argument is either 'method' or list with entry 'method'!");
  if (!method_desc)
    rerror ("ent_jvm_CallMethod: missing descriptor for the method!");

  if (!strlen(meth0))
    rerror ("ent_jvm_CallMethod: string 'method' cannot be of zero-length!");

  // process 'method'
  // two cases:
  //  method="class1.class2.class3.....method
  // or
  //  method="class1/class2/..../method
  //  TODO: method="class1/class2...classK:field:method"
  method = cpstr( meth0 );
  for (i=strlen(method)-2; i>0;i--)
  {
    if(method[i]=='.' || method[i]=='/')
    {
      // convert to JNI convention
      if (method[i]=='.')
        method[i]='/';
      // last '/' delineates the name of the method
      if (!hn)
      {
        hn = i;
        method[hn] = '\0';
        method_name  = &method[hn+1];
        method_class =  method;
      }
    }
  }
  if (!hn)
    rerror ("ent_jvm_CallMethod: string 'method' cannot be of zero-length!");

  // process method.descriptor if it is provided
  // and figure out the return type of method
  for (i=strlen(method_desc)-2; i>0;i--)
  {
    if(method_desc[i]==')')
    {
      desc_rval = &method_desc[i+1];
      break;
    }
  }
  if (!desc_rval)
    rerror ("ent_jvm_CallMethod: method descriptor misses the return type!");

  rent = ent_Create ();

  // Do we need this?
//   res = (*jvm)->AttachCurrentThread(jvm, (void **)&env, NULL);
//   if (res < 0)
//   {
//     fprintf(stderr, "ent_jvm_CallMethod: Thread attach failed!\n");
//     goto jvm_CallMethod_out;
//   }

  // is there such a class?
  jCls = (*env)->FindClass(env, method_class);
  if (jCls == NULL)
  {
    fprintf(stderr, "ent_jvm_CallMethod: Can't find %s class!\n", method_class);

    // error: empty return
    ent_data (rent) = mdr_Create(0,0);
    ent_type (rent) = MATRIX_DENSE_REAL;
    goto jvm_CallMethod_out;
  }

  // construct the descriptor for the method from the list of arguments
  // is there such a method?
  // TODO?  we could construct a descriptor from input data but user would still have to provide
  //        return type of the method
  // in JNI the descriptor does everything: input parameters and their order, and the return type
  jMid = (*env)->GetMethodID(env, jCls, method_name, method_desc);
  if (jMid == NULL)
  {
    jMid = (*env)->GetStaticMethodID(env, jCls, method_name, method_desc);
    ist_meth = 1;
    if (jMid == NULL)
    {
      fprintf(stderr, "ent_jvm_CallMethod: Can't find %s/%s method\n", method_class, method_name);
      fprintf(stderr, "ent_jvm_CallMethod: with descriptor %s!\n", method_desc);
      // error: empty return
      ent_data (rent) = mdr_Create(0,0);
      ent_type (rent) = MATRIX_DENSE_REAL;
      goto jvm_CallMethod_out;
    }
  }

  // if we are here, then the class/method has been found and is ready to be used.
  // we process the arguments and attach them to the jObject array that will pass
  // them to the method being called
  if (nargs>1)
  {
    jArgs = (jvalue *) GC_malloc((nargs-1)*sizeof(jvalue));

    // determine number of arguments: for now accept only strings and int/reals
    for (i=1;i<nargs;i++)
    {
      es = bltin_get_ent (args[i]);

      if (ent_type(es) == MATRIX_DENSE_STRING)
      {
        s  = ent_data(es);
        ns = s->nrow * s->ncol;
        jObj = (*env)->NewObjectArray(env, ns, (*env)->FindClass(env, "java/lang/String"), NULL);
        if(jObj == NULL)
        {
          fprintf(stderr, "ent_jvm_CallMethod: Cannot create object from argument %i\n", i+1);
          // error: empty return
          ent_data (rent) = mdr_Create(0,0);
          ent_type (rent) = MATRIX_DENSE_REAL;
          goto jvm_CallMethod_out;
        }
        for (j=0; j<ns;j++)
        {
          jstring js = (*env)->NewStringUTF(env, MdsV0(s,j));
          (*env)->SetObjectArrayElement(env, jObj, j, js);
//           (*env)->DeleteLocalRef(env, js);
        }
        // clean up 'js' ?
        //
      }
      else if (ent_type(es) == MATRIX_DENSE_REAL)
      {
        r  = ent_data(es);
        if (r->nrow == 1 || r->ncol==1)
        {
          // it is a vector
          ns = r->nrow * r->ncol;
          switch(r->type)
          {
            case RLAB_TYPE_DOUBLE:
              jObj = (*env)->NewDoubleArray(env, ns);
              if(jObj == NULL)
              {
                fprintf(stderr, "ent_jvm_CallMethod: Cannot create object from argument %i\n", i+1);
                // error: empty return
                ent_Clean(es);
                ent_data (rent) = mdr_Create(0,0);
                ent_type (rent) = MATRIX_DENSE_REAL;
                goto jvm_CallMethod_out;
              }
              (*env)->SetDoubleArrayRegion(env, jObj, 0, ns, (jdouble *)MDRPTR(r));
              break;

            case RLAB_TYPE_INT32:
              jObj = (*env)->NewIntArray(env, ns);
              if(jObj == NULL)
              {
                fprintf(stderr, "ent_jvm_CallMethod: Cannot create object from argument %i\n", i+1);
                // error: empty return
                ent_Clean(es);
                ent_data (rent) = mdr_Create(0,0);
                ent_type (rent) = MATRIX_DENSE_REAL;
                goto jvm_CallMethod_out;
              }
              (*env)->SetIntArrayRegion(env, jObj, 0, ns, (jint *) MDIPTR(r));
              break;

            default:
              break;
          }
        }
        else
        {
          // it is a matrix
        }

      }
      else
        printf("ent_jvm_CallMethod: the data type %s is not accepted\n", etd(es));

      jArgs[i-1].l = jObj;
      ent_Clean(es);

    } // for (i=1;i<nargs;i++)

  } // if (nargs>1)

  // call the method
  switch(desc_rval[0])
  {
    case 'V':
      // void method
      if(ist_meth)
        (*env)->CallStaticVoidMethodA(env, jCls, jMid, jArgs);
      else
        (*env)->CallVoidMethodA(env, jCls, jMid, jArgs);

      iretval = 0;
      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'D':
      // double method
      if(ist_meth)
        retval = (*env)->CallStaticDoubleMethodA(env, jCls, jMid, jArgs);
      else
        retval = (*env)->CallDoubleMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdr_CreateScalar(retval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'F':
      // float method
      if(ist_meth)
        retval = (*env)->CallStaticFloatMethodA(env, jCls, jMid, jArgs);
      else
        retval = (*env)->CallFloatMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdr_CreateScalar(retval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'J':
      // long method
      if(ist_meth)
        iretval = (*env)->CallStaticLongMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallLongMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'I':
      // int method
      if(ist_meth)
        iretval = (*env)->CallStaticIntMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallIntMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'S':
      // short method
      if(ist_meth)
        iretval = (*env)->CallStaticShortMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallShortMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'C':
      // char method
      if(ist_meth)
        iretval = (*env)->CallStaticCharMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallCharMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'B':
      // byte
      if(ist_meth)
        iretval = (*env)->CallStaticByteMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallByteMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'Z':
      // boolean
      if(ist_meth)
        iretval = (*env)->CallStaticBooleanMethodA(env, jCls, jMid, jArgs);
      else
        iretval = (*env)->CallBooleanMethodA(env, jCls, jMid, jArgs);

      ent_data (rent) = mdi_CreateScalar(iretval);
      ent_type (rent) = MATRIX_DENSE_REAL;
      break;

    case 'L':
      // string
      if (!strcmp(desc_rval,"Ljava/lang/String;"))
      {
        // string value
        if(ist_meth)
          jRetarray = (*env)->CallStaticObjectMethodA(env, jCls, jMid, jArgs);
        else
          jRetarray = (*env)->CallObjectMethodA(env, jCls, jMid, jArgs);

        int k=(*env)->GetStringLength(env,jRetarray);
        char *c = GC_malloc((k+1)*sizeof(char));
        for (i=0;i<k;i++)
          (*env)->GetStringRegion(env, jRetarray, i, 1, (jchar *)&c[i]);
        c[k]='\0';
//      why the heck this does not work? it returns only the first character of the string
//         (*env)->GetStringRegion(env, jRetarray, 0, k, (jchar *) &c[0]);
        sval = mds_CreateScalar( c );
        ent_data (rent) = sval;
        ent_type (rent) = MATRIX_DENSE_STRING;
      }
      else
      {
        fprintf(stderr, "ent_jvm_CallMethod: Don't know what to do about method return value %s.\n", desc_rval);
        // error: empty return
        ent_data (rent) = mdr_Create(0,0);
        ent_type (rent) = MATRIX_DENSE_REAL;
        goto jvm_CallMethod_out;
      }
      break;

    case '[':
      //
      // we assume that desc_rval[1] exists
      //
      if (!desc_rval[1])
      {
        fprintf(stderr, "ent_jvm_CallMethod: Incorrect descriptor of return value %s. Array of what?\n", desc_rval);
        // error: empty return
        ent_data (rent) = mdr_Create(0,0);
        ent_type (rent) = MATRIX_DENSE_REAL;
        goto jvm_CallMethod_out;

      }

      // array of objects
      if(ist_meth)
        jRetarray = (*env)->CallStaticObjectMethodA(env, jCls, jMid, jArgs);
      else
        jRetarray = (*env)->CallObjectMethodA(env, jCls, jMid, jArgs);

      // figure out the size
      nr = (int)(*env)->GetArrayLength(env, jRetarray);
      if (!nr)
      {
        rval = mdr_Create(nr, 1);
        ent_data (rent) = rval;
        ent_type (rent) = MATRIX_DENSE_REAL;
        break;
      }

      // from the descriptor figure out the return values
      switch(desc_rval[1])
      {
        case 'D':
        case 'F':
        case 'J':
          rval = mdr_Create(nr, 1);
          switch (desc_rval[1])
          {
            case 'D':
              (*env)->GetDoubleArrayRegion(env, jRetarray, 0, nr, MDRPTR(rval));
              break;
            case 'F':
              (*env)->GetFloatArrayRegion(env, jRetarray, 0, nr, (jfloat *)MDRPTR(rval));
              // size(float)=4, while size(double)=8
              // so we need to expand float array to double
              jfloat *jf = (jfloat *)(MDRPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdrV0(rval,i) = jf[i];
              break;
            case 'J':
              (*env)->GetLongArrayRegion(env, jRetarray, 0, nr, (jlong *)MDRPTR(rval));
              jlong *jl = (jlong *)(MDRPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdrV0(rval,i) = jl[i];
              break;
          }
          ent_data (rent) = rval;
          ent_type (rent) = MATRIX_DENSE_REAL;
          break;

        case 'S':
        case 'I':
        case 'C':
        case 'B':
        case 'Z':
          rval = mdi_Create(nr, 1);
          switch (desc_rval[1])
          {
            case 'S':
              (*env)->GetShortArrayRegion(env, jRetarray, 0, nr, (jshort *)MDIPTR(rval));
              jshort *js = (jshort *)(MDIPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdiV0(rval,i) = js[i];
              break;
            case 'I':
              (*env)->GetIntArrayRegion(env, jRetarray, 0, nr, (jint *)MDIPTR(rval));
              break;
            case 'C':
              (*env)->GetCharArrayRegion(env, jRetarray, 0, nr, (jchar*)MDIPTR(rval));
              jchar *jc = (jchar *)(MDIPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdiV0(rval,i) = jc[i];
              break;
            case 'B':
              (*env)->GetByteArrayRegion(env, jRetarray, 0, nr, (jbyte*)MDIPTR(rval));
              jbyte *jb = (jbyte *)(MDIPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdiV0(rval,i) = jb[i];
              break;
            case 'Z':
              (*env)->GetBooleanArrayRegion(env, jRetarray, 0, nr, (jboolean*)MDIPTR(rval));
              jboolean *jz = (jboolean *)(MDIPTR(rval));
              for(i=nr-1;i>=0;i--)
                MdiV0(rval,i) = jz[i];
              break;
          }
          ent_data (rent) = rval;
          ent_type (rent) = MATRIX_DENSE_REAL;
          break;

        case 'L':
          // string
          if (!strcmp(desc_rval,"[Ljava/lang/String;"))
          {
            sval = mds_Create(nr, 1);
            for (i=0; i<nr; i++)
            {
              jobject jo = (*env)->GetObjectArrayElement(env, jRetarray, i);
              int j,k=(*env)->GetStringLength(env,jo);
              char *c = GC_malloc((k+1)*sizeof(char));
              for (j=0;j<k;j++)
                (*env)->GetStringRegion(env, jo, j, 1, (jchar *)&c[j]);
              c[k]='\0';
              MdsV0(sval,i) = c;
            }
            ent_data (rent) = sval;
            ent_type (rent) = MATRIX_DENSE_STRING;
            break;
          }
          else
          {
            fprintf(stderr, "ent_jvm_CallMethod: Don't know what to do about method return value %s.\n", desc_rval);
            // error: empty return
            ent_data (rent) = mdr_Create(0,0);
            ent_type (rent) = MATRIX_DENSE_REAL;
            goto jvm_CallMethod_out;
          }
          break;

      } // switch(desc_rval[1])

      (*env)->DeleteLocalRef(env, jRetarray);
      break;

  } // switch(desc_rval[0])

jvm_CallMethod_out:

  // clean_up
  if (jArgs)
    GC_free (jArgs);
  if (method)
    GC_free (method);
  ent_Clean (e1);

  return (rent);
}

#endif  /* HAVE_JAVA2 */
