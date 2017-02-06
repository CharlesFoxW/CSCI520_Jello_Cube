/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <string>

#define UNDEFORMED_STRUCT_LEN    0.142857

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {

                point netForce = {};
                point normal = {};
                point collisionForce = {};
                int roundedX, roundedY, roundedZ;


                //structural springs: (6)
                if (i+1 >= 0 && i+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i+1][j][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j][k],
                                             jello->v[i][j][k], jello->v[i+1][j][k]), netForce);
                }

                if (j+1 >= 0 && j+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i][j+1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j+1][k],
                                             jello->v[i][j][k], jello->v[i][j+1][k]), netForce);
                }
                if (k+1 >= 0 && k+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i][j][k+1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j][k+1],
                                             jello->v[i][j][k], jello->v[i][j][k+1]), netForce);
                }
                if (i-1 >= 0 && i-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i-1][j][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j][k],
                                             jello->v[i][j][k], jello->v[i-1][j][k]), netForce);
                }
                if (j-1 >= 0 && j-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i][j-1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j-1][k],
                                             jello->v[i][j][k], jello->v[i][j-1][k]), netForce);
                }
                if (k-1 >= 0 && k-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN,
                                             jello->p[i][j][k], jello->p[i][j][k-1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j][k-1],
                                             jello->v[i][j][k], jello->v[i][j][k-1]), netForce);
                }

                //shear springs: (20)
                //12 plain diagonals:
                if (i+1 >= 0 && i+1 <=7 && j+1 >= 0 && j+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i+1][j+1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j+1][k],
                                             jello->v[i][j][k], jello->v[i+1][j+1][k]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j+1 >= 0 && j+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i - 1][j + 1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j+1][k],
                                             jello->v[i][j][k], jello->v[i-1][j+1][k]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j-1 >= 0 && j-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i - 1][j - 1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j-1][k],
                                             jello->v[i][j][k], jello->v[i-1][j-1][k]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && j-1 >= 0 && j-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i + 1][j - 1][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j-1][k],
                                             jello->v[i][j][k], jello->v[i+1][j-1][k]), netForce);
                }
                if (j+1 >= 0 && j+1 <=7 && k+1 >= 0 && k+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i][j + 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j+1][k+1],
                                             jello->v[i][j][k], jello->v[i][j+1][k+1]), netForce);
                }
                if (j-1 >= 0 && j-1 <=7 && k+1 >= 0 && k+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i][j - 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j-1][k+1],
                                             jello->v[i][j][k], jello->v[i][j-1][k+1]), netForce);
                }
                if (j-1 >= 0 && j-1 <=7 && k-1 >= 0 && k-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i][j - 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j-1][k-1],
                                             jello->v[i][j][k], jello->v[i][j-1][k-1]), netForce);
                }
                if (j+1 >= 0 && j+1 <=7 && k-1 >= 0 && k-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i][j + 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j+1][k-1],
                                             jello->v[i][j][k], jello->v[i][j+1][k-1]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && k+1 >= 0 && k+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i + 1][j][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j][k+1],
                                             jello->v[i][j][k], jello->v[i+1][j][k+1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && k+1 >= 0 && k+1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i - 1][j][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j][k+1],
                                             jello->v[i][j][k], jello->v[i-1][j][k+1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && k-1 >= 0 && k-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i - 1][j][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j][k-1],
                                             jello->v[i][j][k], jello->v[i-1][j][k-1]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && k-1 >= 0 && k-1 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(2),
                                             jello->p[i][j][k], jello->p[i + 1][j][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j][k-1],
                                             jello->v[i][j][k], jello->v[i+1][j][k-1]), netForce);
                }

                //8 cube diagonals:
                if (i+1 >= 0 && i+1 <=7 && j+1 >= 0 && j+1 <=7 && k+1 >= 0 && k+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j+1][k+1],
                                             jello->v[i][j][k], jello->v[i+1][j+1][k+1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j+1 >= 0 && j+1 <=7 && k+1 >= 0 && k+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j+1][k+1],
                                             jello->v[i][j][k], jello->v[i-1][j+1][k+1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j-1 >= 0 && j-1 <=7 && k+1 >= 0 && k+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j-1][k+1],
                                             jello->v[i][j][k], jello->v[i-1][j-1][k+1]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && j-1 >= 0 && j-1 <=7 && k+1 >= 0 && k+1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j-1][k+1],
                                             jello->v[i][j][k], jello->v[i+1][j-1][k+1]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && j+1 >= 0 && j+1 <=7 && k-1 >= 0 && k-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j+1][k-1],
                                             jello->v[i][j][k], jello->v[i+1][j+1][k-1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j+1 >= 0 && j+1 <=7 && k-1 >= 0 && k-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j+1][k-1],
                                             jello->v[i][j][k], jello->v[i-1][j+1][k-1]), netForce);
                }
                if (i-1 >= 0 && i-1 <=7 && j-1 >= 0 && j-1 <=7 && k-1 >= 0 && k-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j-1][k-1],
                                             jello->v[i][j][k], jello->v[i-1][j-1][k-1]), netForce);
                }
                if (i+1 >= 0 && i+1 <=7 && j-1 >= 0 && j-1 <=7 && k-1 >= 0 && k-1 <= 7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * sqrt(3),
                                             jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+1][j-1][k-1],
                                             jello->v[i][j][k], jello->v[i+1][j-1][k-1]), netForce);
                }

                //bend springs: (6)
                if (i+2 >= 0 && i+2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i + 2][j][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i+2][j][k],
                                             jello->v[i][j][k], jello->v[i+2][j][k]), netForce);
                }
                if (j+2 >= 0 && j+2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i][j + 2][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j+2][k],
                                             jello->v[i][j][k], jello->v[i][j+2][k]), netForce);
                }
                if (k+2 >= 0 && k+2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i][j][k + 2]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j][k+2],
                                             jello->v[i][j][k], jello->v[i][j][k+2]), netForce);
                }
                if (i-2 >= 0 && i-2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i - 2][j][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i-1][j][k],
                                             jello->v[i][j][k], jello->v[i-1][j][k]), netForce);
                }
                if (j-2 >= 0 && j-2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i][j - 2][k]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j-2][k],
                                             jello->v[i][j][k], jello->v[i][j-2][k]), netForce);
                }
                if (k-2 >= 0 && k-2 <=7) {
                    pSUM(netForce, hookForce(jello->kElastic, UNDEFORMED_STRUCT_LEN * 2,
                                             jello->p[i][j][k], jello->p[i][j][k - 2]), netForce);
                    pSUM(netForce, dampForce(jello->dElastic, jello->p[i][j][k], jello->p[i][j][k-2],
                                             jello->v[i][j][k], jello->v[i][j][k-2]), netForce);
                }

                //collision detection:

                if (jello->p[i][j][k].x < -2.0) {
                    normal.x = 1;
                    normal.y = 0;
                    normal.z = 0;
                    collisionForce.x += jello->kCollision * fabs(jello->p[i][j][k].x + 2.0) * normal.x;
                    collisionForce.x += -1.0 * jello->dCollision * jello->v[i][j][k].x;
                }
                if (jello->p[i][j][k].x > 2.0) {
                    normal.x = -1;
                    normal.y = 0;
                    normal.z = 0;
                    collisionForce.x += jello->kCollision * fabs(jello->p[i][j][k].x - 2.0) * normal.x;
                    collisionForce.x += -1.0 * jello->dCollision * jello->v[i][j][k].x;
                }
                if (jello->p[i][j][k].y < -2.0) {
                    normal.x = 0;
                    normal.y = 1;
                    normal.z = 0;
                    collisionForce.y += jello->kCollision * fabs(jello->p[i][j][k].y + 2.0) * normal.y;
                    collisionForce.y += -1.0 * jello->dCollision * jello->v[i][j][k].y;
                }
                if (jello->p[i][j][k].y > 2.0) {
                    normal.x = 0;
                    normal.y = -1;
                    normal.z = 0;
                    collisionForce.y += jello->kCollision * fabs(jello->p[i][j][k].y - 2.0) * normal.y;
                    collisionForce.y += -1.0 * jello->dCollision * jello->v[i][j][k].y;
                }
                if (jello->p[i][j][k].z < -2.0) {
                    normal.x = 0;
                    normal.y = 0;
                    normal.z = 1;
                    collisionForce.z += jello->kCollision * fabs(jello->p[i][j][k].z + 2.0) * normal.z;
                    collisionForce.z += -1.0 * jello->dCollision * jello->v[i][j][k].z;
                }
                if (jello->p[i][j][k].z > 2.0) {
                    normal.x = 0;
                    normal.y = 0;
                    normal.z = -1;
                    collisionForce.z += jello->kCollision * fabs(jello->p[i][j][k].z - 2.0) * normal.z;
                    collisionForce.z += -1.0 * jello->dCollision * jello->v[i][j][k].z;
                }

                pSUM(netForce, collisionForce, netForce);

                //external force field:
                if (jello->resolution != 0) {

                    int xPos, yPos, zPos;
                    double gridSize = 4.0 / jello->resolution;

                    if (jello->p[i][j][k].x <= -2.0)
                        xPos = 0;
                    else if (jello->p[i][j][k].x >= 2.0)
                        xPos = (int) (4.0 / gridSize) - 1;
                    else
                        xPos = (int) ((jello->p[i][j][k].x + 2.0) / gridSize);

                    if (jello->p[i][j][k].y <= -2.0)
                        yPos = 0;
                    else if (jello->p[i][j][k].y >= 2.0)
                        yPos = (int) (4.0 / gridSize) - 1;
                    else
                        yPos = (int) ((jello->p[i][j][k].y + 2.0) / gridSize);

                    if (jello->p[i][j][k].z <= -2.0)
                        zPos = 0;
                    else if (jello->p[i][j][k].z >= 2.0)
                        zPos = (int) (4.0 / gridSize) - 1;
                    else
                        zPos = (int) ((jello->p[i][j][k].z + 2.0) / gridSize);

                    pSUM(netForce, jello->forceField[xPos * jello->resolution * jello->resolution
                                                     + yPos * jello->resolution +zPos], netForce);

                }


                roundedX = (int) (netForce.x * 100.0);
                roundedY = (int) (netForce.y * 100.0);
                roundedZ = (int) (netForce.z * 100.0);
                netForce.x = (double) roundedX / 100.0;
                netForce.y = (double) roundedY / 100.0;
                netForce.z = (double) roundedZ / 100.0;

                //printf("x = %f\n", netForce.x);
                //printf("y = %f\n", netForce.y);
                //printf("z = %f\n", netForce.z);

                pMULTIPLY(netForce, 1.0 / jello->mass, a[i][j][k]);
            }
        }
    }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);



  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}

double distance(struct point a, struct point b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

point hookForce(double k, double rest, struct point a, struct point b) {
    point dest = {};
    dest.x = -1.0 * k * (distance(a, b) - rest) * (a.x - b.x) / distance(a, b);
    dest.y = -1.0 * k * (distance(a, b) - rest) * (a.y - b.y) / distance(a, b);
    dest.z = -1.0 * k * (distance(a, b) - rest) * (a.z - b.z) / distance(a, b);
    return dest;
}

point dampForce(double d, struct point a, struct point b, struct point va, struct point vb) {
    point dest = {};
    point vDiff = {};
    pDIFFERENCE(va, vb, vDiff);
    double coeff = (vDiff.x * (a.x - b.x) + vDiff.y * (a.y - b.y) + vDiff.z * (a.z - b.z)) / distance(a, b);
    dest.x = -1.0 * d * coeff * (a.x - b.x) / distance(a, b);
    dest.y = -1.0 * d * coeff * (a.y - b.y) / distance(a, b);
    dest.z = -1.0 * d * coeff * (a.z - b.z) / distance(a, b);
    return dest;
}

