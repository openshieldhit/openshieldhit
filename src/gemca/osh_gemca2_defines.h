#ifndef _OSH_GEMCA2_DEFINES
#define _OSH_GEMCA2_DEFINES

/* maximum length of a user given body name, must be sufficiently less than OSH_MAX_LINE_LENGTH:
    NBODY * OSH_GEMCA_BODY_NAME_MAXLEN < OSH_MAX_LINE_LENGTH
 */
#define OSH_GEMCA_BODY_NAME_MAXLEN 0xFF

/* list number of arguments needed to describe each body */
#define OSH_GEMCA_NARGS_MAX 24 /* maximum number of arguments */

#define OSH_GEMCA_NARGS_SPH 4
#define OSH_GEMCA_NARGS_WED 12
#define OSH_GEMCA_NARGS_ARB 24
#define OSH_GEMCA_NARGS_BOX 12
#define OSH_GEMCA_NARGS_VOX 6
#define OSH_GEMCA_NARGS_RPP 6
#define OSH_GEMCA_NARGS_RCC 7
#define OSH_GEMCA_NARGS_REC 12
#define OSH_GEMCA_NARGS_TRC 8
#define OSH_GEMCA_NARGS_ELL 12
#define OSH_GEMCA_NARGS_MOV 4 // TODO
#define OSH_GEMCA_NARGS_ROT 3 // TODO

/* body IDs */
#define OSH_GEMCA_BODY_NONE 0
#define OSH_GEMCA_BODY_SPH 1
#define OSH_GEMCA_BODY_WED 2
#define OSH_GEMCA_BODY_ARB 3

#define OSH_GEMCA_BODY_BOX 4
#define OSH_GEMCA_BODY_VOX 5
#define OSH_GEMCA_BODY_RPP 6

#define OSH_GEMCA_BODY_RCC 7
#define OSH_GEMCA_BODY_REC 8
#define OSH_GEMCA_BODY_TRC 9
#define OSH_GEMCA_BODY_ELL 10

#define OSH_GEMCA_BODY_YZP 11
#define OSH_GEMCA_BODY_XZP 12
#define OSH_GEMCA_BODY_XYP 13
#define OSH_GEMCA_BODY_PLA 14

#define OSH_GEMCA_BODY_ROT 100 // TODO
#define OSH_GEMCA_BODY_CPY 101 // TODO
#define OSH_GEMCA_BODY_MOV 102 // TODO

#define OSH_GEMCA_BODY_LOGIC_NONE 0
#define OSH_GEMCA_BODY_LOGIC_INCL 1 /* eg. "+1" */
#define OSH_GEMCA_BODY_LOGIC_EXCL 2 /* eg. "-1" */
#define OSH_GEMCA_BODY_LOGIC_OR 3   /* eg. "OR1"*/

/* surface types */
#define OSH_GEMCA_SURF_NONE 0
#define OSH_GEMCA_SURF_SPHERE 1    /* 3D surface */
#define OSH_GEMCA_SURF_ELLIPSOID 2 /* 3D surface */
#define OSH_GEMCA_SURF_CYLZ 3      /* inifinte cylinder along Z, defined as a circle in 2D x-y plane */
#define OSH_GEMCA_SURF_ELLZ 4      /* infinite elliptical cylinder along Z, defined as an elliplse in 2D x-y plane */
#define OSH_GEMCA_SURF_CONE 5      /* inifinite cone surface */
#define OSH_GEMCA_SURF_PLANEX 6    /* a YZ plane, defined by Ax - D = 0 */
#define OSH_GEMCA_SURF_PLANEY 7    /* a XZ plane, defined by By - D = 0 */
#define OSH_GEMCA_SURF_PLANEZ 8    /* a XY plane, defined by Cz - D = 0 */
#define OSH_GEMCA_SURF_PLANE 9     /* a plane arbitrary oriented and located in space, Ax + Bx + Cz - D = 0 */

#define OSH_GEMCA_SURF_POSITIVE 1 /* s->sign in positve direction, i.e. normal vector points into the body */
#define OSH_GEMCA_SURF_NEGATIVE 0 /* s->sign in negtive direction, i.e. normal vector points outside the body */

#endif /* _OSH_GEMCA2_DEFINES */
