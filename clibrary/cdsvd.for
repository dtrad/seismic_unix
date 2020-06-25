
      SUBROUTINE CDSVD(CA,S,U,V,MDIM,NDIM,M,N,NS)                               
      COMPLEX*16 CA(MDIM,N),S(N),U(MDIM,N),V(NDIM,N)                            
      EXTERNAL CDSVDW                                                           
      NS=M                                                                      
      IF(N.LT.M)NS=N                                                            
      CALL GSPACE(XR,M*N*8)                                                     
      CALL GSPACE(XI,M*N*8)                                                     
      CALL GSPACE(UR,M*NS*8)                                                    
      CALL GSPACE(UI,M*NS*8)                                                    
      CALL GSPACE(VR,N*N*8)                                                     
      CALL GSPACE(VI,N*N*8)                                                     
      CALL GSPACE(SR,(NS+1)*8)                                                  
      CALL GSPACE(SI,(NS+1)*8)                                                  
      CALL GSPACE(ER,N*8)                                                       
      CALL GSPACE(EI,N*8)                                                       
      CALL GSPACE(WORKR,M*8)                                                    
      CALL GSPACE(WORKI,M*8)                                                    
      CALL CALLER(CDSVDW,IPTR(CA),IPTR(S),IPTR(U),IPTR(V),                      
     *IPTR(MDIM),IPTR(NDIM),IPTR(M),IPTR(N),IPTR(NS),                           
     *XR,XI,UR,UI,VR,VI,SR,SI,ER,EI,WORKR,WORKI)                                
      CALL FSPACE(XR)                                                           
      CALL FSPACE(XI)                                                           
      CALL FSPACE(UR)                                                           
      CALL FSPACE(UI)                                                           
      CALL FSPACE(VR)                                                           
      CALL FSPACE(VI)                                                           
      CALL FSPACE(SR)                                                           
      CALL FSPACE(SI)                                                           
      CALL FSPACE(ER)                                                           
      CALL FSPACE(EI)                                                           
      CALL FSPACE(WORKR)                                                        
      CALL FSPACE(WORKI)                                                        
      RETURN                                                                    
      END                                                                       
      SUBROUTINE CDSVDW(CA,S,U,V,MDIM,NDIM,M,N,NS,                              
     *XR,XI,UR,UI,VR,VI,SR,SI,ER,EI,WORKR,WORKI)                                
      COMPLEX*16 CA(MDIM,N),S(N),U(MDIM,N),V(NDIM,N)                            
      REAL*8 XR(M,N),XI(M,N),SR(1),SI(1),ER(N),EI(N)                            
      REAL*8 UR(M,N),UI(M,N),VR(N,N),VI(N,N),WORKR(M),WORKI(M)                  
      INTEGER LDX,P,LDU,LDV,JOB,INFO                                            
      REAL*8 DREAL,DIMAG                                                        
      LDX=M                                                                     
      LDU=M                                                                     
      LDV=N                                                                     
      JOB=21                                                                    
      DO 100 I=1,M                                                              
      DO 100 J=1,N                                                              
      XR(I,J)=DREAL(CA(I,J))                                                    
      XI(I,J)=DIMAG(CA(I,J))                                                    
  100 CONTINUE                                                                  
      CALL WSVDC(XR,XI,LDX,M,N,SR,SI,ER,EI,UR,UI,LDU,VR,VI,LDV,                 
     *WORKR,WORKI,JOB,INFO)                                                     
      NS=M                                                                      
      IF(N.LT.M)NS=N                                                            
      DO 200 J=1,NS                                                             
      S(J)=DCMPLX(SR(J),SI(J))                                                  
      DO 200 I=1,M                                                              
  200 U(I,J)=DCMPLX(UR(I,J),UI(I,J))                                            
      DO 300 J=1,N                                                              
      DO 300 I=1,N                                                              
  300 V(I,J)=DCMPLX(VR(I,J),VI(I,J))                                            
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WSVDC(XR,XI,LDX,N,P,SR,SI,ER,EI,UR,UI,LDU,VR,VI,LDV,           
     *                 WORKR,WORKI,JOB,INFO)                                    
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO                                          
      DOUBLE PRECISION XR(LDX,1),XI(LDX,1),SR(1),SI(1),ER(1),EI(1),             
     *                 UR(LDU,1),UI(LDU,1),VR(LDV,1),VI(LDV,1),                 
     *                 WORKR(1),WORKI(1)                                        
C                                                                               
C                                                                               
C     WSVDC IS A SUBROUTINE TO REDUCE A DOUBLE-COMPLEX NXP MATRIX X BY          
C     UNITARY TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE                    
C     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE                 
C     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,                 
C     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.                          
C                                                                               
C     ON ENTRY                                                                  
C                                                                               
C         X         DOUBLE-COMPLEX(LDX,P), WHERE LDX.GE.N.                      
C                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE                  
C                   DECOMPOSITION IS TO BE COMPUTED.  X IS                      
C                   DESTROYED BY WSVDC.                                         
C                                                                               
C         LDX       INTEGER.                                                    
C                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.                
C                                                                               
C         N         INTEGER.                                                    
C                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.                 
C                                                                               
C         P         INTEGER.                                                    
C                   P IS THE NUMBER OF ROWS OF THE MATRIX X.                    
C                                                                               
C         LDU       INTEGER.                                                    
C                   LDU IS THE LEADING DIMENSION OF THE ARRAY U                 
C                   (SEE BELOW).                                                
C                                                                               
C         LDV       INTEGER.                                                    
C                   LDV IS THE LEADING DIMENSION OF THE ARRAY V                 
C                   (SEE BELOW).                                                
C                                                                               
C         WORK      DOUBLE-COMPLEX(N).                                          
C                   WORK IS A SCRATCH ARRAY.                                    
C                                                                               
C         JOB       INTEGER.                                                    
C                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR                
C                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB                   
C                   WITH THE FOLLOWING MEANING                                  
C                                                                               
C                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR             
C                                  VECTORS.                                     
C                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS           
C                                  IN U.                                        
C                        A.GE.2    RETURNS THE FIRST MIN(N,P)                   
C                                  LEFT SINGULAR VECTORS IN U.                  
C                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR            
C                                  VECTORS.                                     
C                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS            
C                                  IN V.                                        
C                                                                               
C     ON RETURN                                                                 
C                                                                               
C         S         DOUBLE-COMPLEX(MM), WHERE MM=MIN(N+1,P).                    
C                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE                 
C                   SINGULAR VALUES OF X ARRANGED IN DESCENDING                 
C                   ORDER OF MAGNITUDE.                                         
C                                                                               
C         E         DOUBLE-COMPLEX(P).                                          
C                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE               
C                   DISCUSSION OF INFO FOR EXCEPTIONS.                          
C                                                                               
C         U         DOUBLE-COMPLEX(LDU,K), WHERE LDU.GE.N.                      
C                   IF JOBA.EQ.1 THEN K.EQ.N,                                   
C                   IF JOBA.EQ.2 THEN K.EQ.MIN(N,P).                            
C                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.            
C                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P                
C                   OR IF JOBA.GT.2, THEN U MAY BE IDENTIFIED WITH X            
C                   IN THE SUBROUTINE CALL.                                     
C                                                                               
C         V         DOUBLE-COMPLEX(LDV,P), WHERE LDV.GE.P.                      
C                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.            
C                   V IS NOT REFERENCED IF JOBB.EQ.0.  IF P.LE.N,               
C                   THEN V MAY BE IDENTIFIED WHTH X IN THE                      
C                   SUBROUTINE CALL.                                            
C                                                                               
C         INFO      INTEGER.                                                    
C                   THE SINGULAR VALUES (AND THEIR CORRESPONDING                
C                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)              
C                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF                     
C                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR                
C                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX              
C                   B = CTRANS(U)*X*V IS THE BIDIAGONAL MATRIX                  
C                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE              
C                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (CTRANS(U)              
C                   IS THE CONJUGATE-TRANSPOSE OF U).  THUS THE                 
C                   SINGULAR VALUES OF X AND B ARE THE SAME.                    
C                                                                               
C     LINPACK. THIS VERSION DATED 07/03/79 .                                    
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.               
C                                                                               
C     WSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.                       
C                                                                               
C     BLAS WAXPY,PYTHAG,WDOTCR,WDOTCI,WSCAL,WSWAP,WNRM2,RROTG                   
C     FORTRAN DABS,DIMAG,DMAX1                                                  
C     FORTRAN MAX0,MIN0,MOD,DSQRT                                               
C                                                                               
C     INTERNAL VARIABLES                                                        
C                                                                               
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,           
     *        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1                                
      DOUBLE PRECISION PYTHAG,WDOTCR,WDOTCI,TR,TI,RR,RI                         
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,WNRM2,SCALE,SHIFT,SL,SM,SN,           
     *                 SMM1,T1,TEST,ZTEST,SMALL,FLOP                            
      LOGICAL WANTU,WANTV                                                       
C                                                                               
      DOUBLE PRECISION ZDUMR,ZDUMI                                              
      DOUBLE PRECISION CABS1                                                    
      CABS1(ZDUMR,ZDUMI) = DABS(ZDUMR) + DABS(ZDUMI)                            
C                                                                               
C     SET THE MAXIMUM NUMBER OF ITERATIONS.                                     
C                                                                               
      MAXIT = 75                                                                
C                                                                               
C     SMALL NUMBER, ROUGHLY MACHINE EPSILON, USED TO AVOID UNDERFLOW            
C                                                                               
      SMALL = 1.D0/2.D0**48                                                     
C                                                                               
C     DETERMINE WHAT IS TO BE COMPUTED.                                         
C                                                                               
      WANTU = .FALSE.                                                           
      WANTV = .FALSE.                                                           
      JOBU = MOD(JOB,100)/10                                                    
      NCU = N                                                                   
      IF (JOBU .GT. 1) NCU = MIN0(N,P)                                          
      IF (JOBU .NE. 0) WANTU = .TRUE.                                           
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.                                    
C                                                                               
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS                
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.                                
C                                                                               
      INFO = 0                                                                  
      NCT = MIN0(N-1,P)                                                         
      NRT = MAX0(0,MIN0(P-2,N))                                                 
      LU = MAX0(NCT,NRT)                                                        
      IF (LU .LT. 1) GO TO 190                                                  
      DO 180 L = 1, LU                                                          
         LP1 = L + 1                                                            
         IF (L .GT. NCT) GO TO 30                                               
C                                                                               
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND                  
C           PLACE THE L-TH DIAGONAL IN S(L).                                    
C                                                                               
            SR(L) = WNRM2(N-L+1,XR(L,L),XI(L,L),1)                              
            SI(L) = 0.0D0                                                       
            IF (CABS1(SR(L),SI(L)) .EQ. 0.0D0) GO TO 20                         
               IF (CABS1(XR(L,L),XI(L,L)) .EQ. 0.0D0) GO TO 10                  
                  CALL WSIGN(SR(L),SI(L),XR(L,L),XI(L,L),SR(L),SI(L))           
   10          CONTINUE                                                         
               CALL WDIV(1.0D0,0.0D0,SR(L),SI(L),TR,TI)                         
               CALL WSCAL(N-L+1,TR,TI,XR(L,L),XI(L,L),1)                        
               XR(L,L) = FLOP(1.0D0 + XR(L,L))                                  
   20       CONTINUE                                                            
            SR(L) = -SR(L)                                                      
            SI(L) = -SI(L)                                                      
   30    CONTINUE                                                               
         IF (P .LT. LP1) GO TO 60                                               
         DO 50 J = LP1, P                                                       
            IF (L .GT. NCT) GO TO 40                                            
            IF (CABS1(SR(L),SI(L)) .EQ. 0.0D0) GO TO 40                         
C                                                                               
C              APPLY THE TRANSFORMATION.                                        
C                                                                               
               TR = -WDOTCR(N-L+1,XR(L,L),XI(L,L),1,XR(L,J),XI(L,J),1)          
               TI = -WDOTCI(N-L+1,XR(L,L),XI(L,L),1,XR(L,J),XI(L,J),1)          
               CALL WDIV(TR,TI,XR(L,L),XI(L,L),TR,TI)                           
               CALL WAXPY(N-L+1,TR,TI,XR(L,L),XI(L,L),1,XR(L,J),                
     *                    XI(L,J),1)                                            
   40       CONTINUE                                                            
C                                                                               
C           PLACE THE L-TH ROW OF X INTO  E FOR THE                             
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.                   
C                                                                               
            ER(J) = XR(L,J)                                                     
            EI(J) = -XI(L,J)                                                    
   50    CONTINUE                                                               
   60    CONTINUE                                                               
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 80                               
C                                                                               
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK                   
C           MULTIPLICATION.                                                     
C                                                                               
            DO 70 I = L, N                                                      
               UR(I,L) = XR(I,L)                                                
               UI(I,L) = XI(I,L)                                                
   70       CONTINUE                                                            
   80    CONTINUE                                                               
         IF (L .GT. NRT) GO TO 170                                              
C                                                                               
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE                   
C           L-TH SUPER-DIAGONAL IN E(L).                                        
C                                                                               
            ER(L) = WNRM2(P-L,ER(LP1),EI(LP1),1)                                
            EI(L) = 0.0D0                                                       
            IF (CABS1(ER(L),EI(L)) .EQ. 0.0D0) GO TO 100                        
               IF (CABS1(ER(LP1),EI(LP1)) .EQ. 0.0D0) GO TO 90                  
                  CALL WSIGN(ER(L),EI(L),ER(LP1),EI(LP1),ER(L),EI(L))           
   90          CONTINUE                                                         
               CALL WDIV(1.0D0,0.0D0,ER(L),EI(L),TR,TI)                         
               CALL WSCAL(P-L,TR,TI,ER(LP1),EI(LP1),1)                          
               ER(LP1) = FLOP(1.0D0 + ER(LP1))                                  
  100       CONTINUE                                                            
            ER(L) = -ER(L)                                                      
            EI(L) = +EI(L)                                                      
            IF (LP1 .GT. N .OR. CABS1(ER(L),EI(L)) .EQ. 0.0D0)                  
     *         GO TO 140                                                        
C                                                                               
C              APPLY THE TRANSFORMATION.                                        
C                                                                               
               DO 110 I = LP1, N                                                
                  WORKR(I) = 0.0D0                                              
                  WORKI(I) = 0.0D0                                              
  110          CONTINUE                                                         
               DO 120 J = LP1, P                                                
                  CALL WAXPY(N-L,ER(J),EI(J),XR(LP1,J),XI(LP1,J),1,             
     *                       WORKR(LP1),WORKI(LP1),1)                           
  120          CONTINUE                                                         
               DO 130 J = LP1, P                                                
                  CALL WDIV(-ER(J),-EI(J),ER(LP1),EI(LP1),TR,TI)                
                  CALL WAXPY(N-L,TR,-TI,WORKR(LP1),WORKI(LP1),1,                
     *                       XR(LP1,J),XI(LP1,J),1)                             
  130          CONTINUE                                                         
  140       CONTINUE                                                            
            IF (.NOT.WANTV) GO TO 160                                           
C                                                                               
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT                     
C              BACK MULTIPLICATION.                                             
C                                                                               
               DO 150 I = LP1, P                                                
                  VR(I,L) = ER(I)                                               
                  VI(I,L) = EI(I)                        LC;∞Ì[∂ßfWóÑL˚(ÄåmC∑‚ï=ã;ÑÉJ»Ô◊6∞öVéÌÍ’∞µhSpﬂWÊ2à∑√s¬≥ƒû‚Õ®˝≥>m´ΩïŸÑU“Hó˝≤I∞–mı£(tnQd3“]«πÏÉC))’!M¥∫)…@úÇ⁄yzÿT§˚¯ı1‚(ÕÁı:G9ÄˇÑıGd:ZYΩ:ÏÑÛ™1?ﬁåüÎçxiµØÑN“∑·ënw7›ÇfΩkıHäÌÄk;ﬁ%¢äΩGJù§◊$"¢’˝≥x3\¨ç^ë–E`Ù?äÍVWW?˝ìbdÎ|òf©|>‘ÂµıeÕ4©'ks˜Ì–#”j.TbaV˝ªm¶ë¢y’%ä}z∞C ø?ÚO∆¶áb-è,‘5ì=Ú ^Ó<Ùæà<†¡`ª&¬¨∏Æ¥ñÏÓ<rÏ\¬%` ø¶ÁÕºù~ﬁ,ÃÓN´ê'æ(Òÿ)Ω@"Å® +Èﬁ}2]$Å⁄‘|…e√îUlíÖvT7á€‰˜û€ú5ÖAOñ>Ÿ bÍz–∞”◊(ÿöä¢Î÷µc–’ﬂ‡ò$_ÉÆ‹Á\≥¸˚y¨Z—∏FÑåîÊΩâèçŸ∂ï‰jíï˚ú5ıb' Ì8›ÄÕ_Q”;ª(5áiŒpö®+±†}øÍzﬂ¯jB|9Gõﬁ±˚∏ÁÀÉU›–   
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.                            
C                                                                               
      M = MIN0(P,N+1)                                                           
      NCTP1 = NCT + 1                                                           
      NRTP1 = NRT + 1                                                           
      IF (NCT .GE. P) GO TO 200                                                 
         SR(NCTP1) = XR(NCTP1,NCTP1)                                            
         SI(NCTP1) = XI(NCTP1,NCTP1)                                            
  200 CONTINUE                                                                  
      IF (N .GE. M) GO TO 210                                                   
         SR(M) = 0.0D0                                                          
         SI(M) = 0.0D0                                                          
  210 CONTINUE                                                                  
      IF (NRTP1 .GE. M) GO TO 220                                               
         ER(NRTP1) = XR(NRTP1,M)                                                
         EI(NRTP1) = XI(NRTP1,M)                                                
  220 CONTINUE                                                                  
      ER(M) = 0.0D0                                                             
      EI(M) = 0.0D0                                                             
C                                                                               
C     IF REQUIRED, GENERATE U.                                                  
C                                                                               
      IF (.NOT.WANTU) GO TO 350                                                 
         IF (NCU .LT. NCTP1) GO TO 250                                          
         DO 240 J = NCTP1, NCU                                                  
            DO 230 I = 1, N                                                     
               UR(I,J) = 0.0D0                                                  
               UI(I,J) = 0.0D0                                                  
  230       CONTINUE                                                            
            UR(J,J) = 1.0D0                                                     
            UI(J,J) = 0.0D0                                                     
  240    CONTINUE                                                               
  250    CONTINUE                                                               
         IF (NCT .LT. 1) GO TO 340                                              
         DO 330 LL = 1, NCT                                                     
            L = NCT - LL + 1                                                    
            IF (CABS1(SR(L),SI(L)) .EQ. 0.0D0) GO TO 300                        
               LP1 = L + 1                                                      
               IF (NCU .LT. LP1) GO TO 270                                      
               DO 260 J = LP1, NCU                                              
                  TR = -WDOTCR(N-L+1,UR(L,L),UI(L,L),1,UR(L,J),                 
     *                         UI(L,J),1)                                       
                  TI = -WDOTCI(N-L+1,UR(L,L),UI(L,L),1,UR(L,J),                 
     *                         UI(L,J),1)                                       
                  CALL WDIV(TR,TI,UR(L,L),UI(L,L),TR,TI)                        
                  CALL WAXPY(N-L+1,TR,TI,UR(L,L),UI(L,L),1,UR(L,J),             
     *                       UI(L,J),1)                                         
  260          CONTINUE                                                         
  270          CONTINUE                                                         
               CALL WRSCAL(N-L+1,-1.0D0,UR(L,L),UI(L,L),1)                      
               UR(L,L) = FLOP(1.0D0 + UR(L,L))                                  
               LM1 = L - 1                                                      
               IF (LM1 .LT. 1) GO TO 290                                        
               DO 280 I = 1, LM1                                                
                  UR(I,L) = 0.0D0                                               
                  UI(I,L) = 0.0D0                                               
  280          CONTINUE                                                         
  290          CONTINUE                                                         
            GO TO 320                                                           
  300       CONTINUE                                                            
               DO 310 I = 1, N                                                  
                  UR(I,L) = 0.0D0                                               
                  UI(I,L) = 0.0D0                                               
  310          CONTINUE                                                         
               UR(L,L) = 1.0D0                                                  
               UI(L,L) = 0.0D0                                                  
  320       CONTINUE                                                            
  330    CONTINUE                                                               
  340    CONTINUE                                                               
  350 CONTINUE                                                                  
C                                                                               
C     IF IT IS REQUIRED, GENERATE V.                                            
C                                                                               
      IF (.NOT.WANTV) GO TO 400                                                 
         DO 390 LL = 1, P                                                       
            L = P - LL + 1                                                      
            LP1 = L + 1                                                         
            IF (L .GT. NRT) GO TO 370                                           
            IF (CABS1(ER(L),EI(L)) .EQ. 0.0D0) GO TO 370                        
               DO 360 J = LP1, P                                                
                  TR = -WDOTCR(P-L,VR(LP1,L),VI(LP1,L),1,VR(LP1,J),             
     *                         VI(LP1,J),1)                                     
                  TI = -WDOTCI(P-L,VR(LP1,L),VI(LP1,L),1,VR(LP1,J),             
     *                         VI(LP1,J),1)                                     
                  CALL WDIV(TR,TI,VR(LP1,L),VI(LP1,L),TR,TI)                    
                  CALL WAXPY(P-L,TR,TI,VR(LP1,L),VI(LP1,L),1,VR(LP1,J),         
     *                       VI(LP1,J),1)                                       
  360          CONTINUE                                                         
  370       CONTINUE                                                            
            DO 380 I = 1, P                                                     
               VR(I,L) = 0.0D0                                                  
               VI(I,L) = 0.0D0                                                  
  380       CONTINUE                                                            
            VR(L,L) = 1.0D0                                                     
            VI(L,L) = 0.0D0                                                     
  390    CONTINUE                                                               
  400 CONTINUE                                                                  
C                                                                               
C     TRANSFORM S AND E SO THAT THEY ARE REAL.                                  
C                                                                               
      DO 420 I = 1, M                                                           
            TR = PYTHAG(SR(I),SI(I))                                            
            IF (TR .EQ. 0.0D0) GO TO 405                                        
            RR = SR(I)/TR                                                       
            RI = SI(I)/TR                                                       
            SR(I) = TR                                                          
            SI(I) = 0.0D0                                                       
            IF (I .LT. M) CALL WDIV(ER(I),EI(I),RR,RI,ER(I),EI(I))              
            IF (WANTU) CALL WSCAL(N,RR,RI,UR(1,I),UI(1,I),1)                    
  405    CONTINUE                                                               
C     ...EXIT                                                                   
         IF (I .EQ. M) GO TO 430                                                
            TR = PYTHAG(ER(I),EI(I))                                            
            IF (TR .EQ. 0.0D0) GO TO 410                                        
            CALL WDIV(TR,0.0D0,ER(I),EI(I),RR,RI)                               
            ER(I) = TR                                                          
            EI(I) = 0.0D0                                                       
            CALL WMUL(SR(I+1),SI(I+1),RR,RI,SR(I+1),SI(I+1))                    
            IF (WANTV) CALL WSCAL(P,RR,RI,VR(1,I+1),VI(1,I+1),1)                
  410    CONTINUE                                                               
  420 CONTINUE                                                                  
  430 CONTINUE                                                                  
C                                                                               
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.                              
C                                                                               
      MM = M                                                                    
      ITER = 0                                                                  
  440 CONTINUE                                                                  
C                                                                               
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.                       
C                                                                               
C     ...EXIT                                                                   
         IF (M .EQ. 0) GO TO 700                                                
C                                                                               
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET                        
C        FLAG AND RETURN.                                                       
C                                                                               
         IF (ITER .LT. MAXIT) GO TO 450                                         
            INFO = M                                                            
C     ......EXIT                                                                
            GO TO 700                                                           
  450    CONTINUE                                                               
C                                                                               
C        THIS SECTION OF THE PROGRAM INSPECTS FOR                               
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON                         
C        COMPLETION THE VARIABLE KASE IS SET AS FOLLOWS.                        
C                                                                               
C           KASE = 1     IF SR(M) AND ER(L-1) ARE NEGLIGIBLE AND L.LT.M         
C           KASE = 2     IF SR(L) IS NEGLIGIBLE AND L.LT.M                      
C           KASE = 3     IF ER(L-1) IS NEGLIGIBLE, L.LT.M, AND                  
C                        SR(L), ..., SR(M) ARE NOT NEGLIGIBLE (QR STEP).        
C           KASE = 4     IF ER(M-1) IS NEGLIGIBLE (CONVERGENCE).                
C                                                                               
         DO 470 LL = 1, M                                                       
            L = M - LL                                                          
C        ...EXIT                                                                
            IF (L .EQ. 0) GO TO 480                                             
            TEST = FLOP(DABS(SR(L)) + DABS(SR(L+1)))                            
            ZTEST = FLOP(TEST + DABS(ER(L))/2.0D0)                              
            IF (SMALL*ZTEST .NE. SMALL*TEST) GO TO 460                          
               ER(L) = 0.0D0                                                    
C        ......EXIT                                                             
               GO TO 480                                                        
  460       CONTINUE                                                            
  470    CONTINUE                                                               
  480    CONTINUE                                                               
         IF (L .NE. M - 1) GO TO 490                                            
            KASE = 4                                                            
         GO TO 560                                                              
  490    CONTINUE                                                               
            LP1 = L + 1                                                         
            MP1 = M + 1                                                         
            DO 510 LLS = LP1, MP1                                               
               LS = M - LLS + LP1                                               
C           ...EXIT                                                             
               IF (LS .EQ. L) GO TO 520                                         
               TEST = 0.0D0                                                     
               IF (LS .NE. M) TEST = FLOP(TEST + DABS(ER(LS)))                  
               IF (LS .NE. L + 1) TEST = FLOP(TEST + DABS(ER(LS-1)))            
               ZTEST = FLOP(TEST + DABS(SR(LS))/2.0D0)                          
               IF (SMALL*ZTEST .NE. SMALL*TEST) GO TO 500                       
                  SR(LS) = 0.0D0                                                
C           ......EXIT                                                          
                  GO TO 520                                                     
  500          CONTINUE                                                         
  510       CONTINUE                                                            
  520       CONTINUE                                                            
            IF (LS .NE. L) GO TO 530                                            
               KASE = 3                                                         
            GO TO 550                                                           
  530       CONTINUE                                                            
            IF (LS .NE. M) GO TO 540                                            
               KASE = 1                                                         
            GO TO 550                                                           
  540       CONTINUE                                                            
               KASE = 2                                                         
               L = LS                                                           
  550       CONTINUE                                                            
  560    CONTINUE                                                               
         L = L + 1                                                              
C                                                                               
C        PERFORM THE TASK INDICATED BY KASE.                                    
C                                                                               
         GO TO (570, 600, 620, 650), KASE                                       
C                                                                               
C        DEFLATE NEGLIGIBLE SR(M).                                              
C                                                                               
  570    CONTINUE                                                               
            MM1 = M - 1                                                         
            F = ER(M-1)                                                         
            ER(M-1) = 0.0D0                                                     
            DO 590 KK = L, MM1                                                  
               K = MM1 - KK + L                                                 
               T1 = SR(K)                                                       
               CALL RROTG(T1,F,CS,SN)                                           
               SR(K) = T1                                                       
               IF (K .EQ. L) GO TO 580                                          
                  F = FLOP(-SN*ER(K-1))                                         
                  ER(K-1) = FLOP(CS*ER(K-1))                                    
  580          CONTINUE                                                         
               IF (WANTV) CALL RROT(P,VR(1,K),1,VR(1,M),1,CS,SN)                
               IF (WANTV) CALL RROT(P,VI(1,K),1,VI(1,M),1,CS,SN)                
  590       CONTINUE                                                            
         GO TO 690                                                              
C                                                                               
C        SPLIT AT NEGLIGIBLE SR(L).                                             
C                                                                               
  600    CONTINUE                                                               
            F = ER(L-1)                                                         
            ER(L-1) = 0.0D0                                                     
            DO 610 K = L, M                                                     
               T1 = SR(K)                                                       
               CALL RROTG(T1,F,CS,SN)                                           
               SR(K) = T1                                                       
               F = FLOP(-SN*ER(K))                                              
               ER(K) = FLOP(CS*ER(K))                                           
               IF (WANTU) CALL RROT(N,UR(1,K),1,UR(1,L-1),1,CS,SN)              
               IF (WANTU) CALL RROT(N,UI(1,K),1,UI(1,L-1),1,CS,SN)              
  610       CONTINUE                                                            
         GO TO 690                                                              
C                                                                               
C        PERFORM ONE QR STEP.                                                   
C                                                                               
  620    CONTINUE                                                               
C                                                                               
C           CALCULATE THE SHIFT.                                                
C                                                                               
            SCALE = DMAX1(DABS(SR(M)),DABS(SR(M-1)),DABS(ER(M-1)),              
     *                    DABS(SR(L)),DABS(ER(L)))                              
            SM = SR(M)/SCALE                                                    
            SMM1 = SR(M-1)/SCALE                                                
            EMM1 = ER(M-1)/SCALE                                                
            SL = SR(L)/SCALE                                                    
            EL = ER(L)/SCALE                                                    
            B = FLOP(((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0)                 
            C = FLOP((SM*EMM1)**2)                                              
            SHIFT = 0.0D0                                                       
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 630                      
               SHIFT = FLOP(DSQRT(B**2+C))                                      
               IF (B .LT. 0.0D0) SHIFT = -SHIFT                                 
               SHIFT = FLOP(C/(B + SHIFT))                                      
  630       CONTINUE                                                            
            F = FLOP((SL + SM)*(SL - SM) - SHIFT)                               
            G = FLOP(SL*EL)                                                     
C                                                                               
C           CHASE ZEROS.                                                        
C                                                                               
            MM1 = M - 1                                                         
            DO 640 K = L, MM1                                                   
               CALL RROTG(F,G,CS,SN)                                            
               IF (K .NE. L) ER(K-1) = F                                        
               F = FLOP(CS*SR(K) + SN*ER(K))                                    
               ER(K) = FLOP(CS*ER(K) - SN*SR(K))                                
               G = FLOP(SN*SR(K+1))                                             
               SR(K+1) = FLOP(CS*SR(K+1))                                       
               IF (WANTV) CALL RROT(P,VR(1,K),1,VR(1,K+1),1,CS,SN)              
               IF (WANTV) CALL RROT(P,VI(1,K),1,VI(1,K+1),1,CS,SN)              
               CALL RROTG(F,G,CS,SN)                                            
               SR(K) = F                                                        
               F = FLOP(CS*ER(K) + SN*SR(K+1))                                  
               SR(K+1) = FLOP(-SN*ER(K) + CS*SR(K+1))                           
               G = FLOP(SN*ER(K+1))                                             
               ER(K+1) = FLOP(CS*ER(K+1))                                       
               IF (WANTU .AND. K .LT. N)                                        
     *            CALL RROT(N,UR(1,K),1,UR(1,K+1),1,CS,SN)                      
               IF (WANTU .AND. K .LT. N)                                        
     *            CALL RROT(N,UI(1,K),1,UI(1,K+1),1,CS,SN)                      
  640       CONTINUE                                                            
            ER(M-1) = F                                                         
            ITER = ITER + 1                                                     
         GO TO 690                                                              
C                                                                               
C        CONVERGENCE                                                            
C                                                                               
  650    CONTINUE                                                               
C                                                                               
C           MAKE THE SINGULAR VALUE  POSITIVE                                   
C                                                                               
            IF (SR(L) .GE. 0.0D0) GO TO 660                                     
               SR(L) = -SR(L)                                                   
             IF (WANTV) CALL WRSCAL(P,-1.0D0,VR(1,L),VI(1,L),1)                 
  660       CONTINUE                                                            
C                                                                               
C           ORDER THE SINGULAR VALUE.                                           
C                                                                               
  670       IF (L .EQ. MM) GO TO 680                                            
C           ...EXIT                                                             
               IF (SR(L) .GE. SR(L+1)) GO TO 680                                
               TR = SR(L)                                                       
               SR(L) = SR(L+1)                                                  
               SR(L+1) = TR                                                     
               IF (WANTV .AND. L .LT. P)                                        
     *            CALL WSWAP(P,VR(1,L),VI(1,L),1,VR(1,L+1),VI(1,L+1),1)         
               IF (WANTU .AND. L .LT. N)                                        
     *            CALL WSWAP(N,UR(1,L),UI(1,L),1,UR(1,L+1),UI(1,L+1),1)         
               L = L + 1                                                        
            GO TO 670                                                           
  680       CONTINUE                                                            
            ITER = 0                                                            
            M = M - 1                                                           
  690    CONTINUE                                                               
      GO TO 440                                                                 
  700 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WAXPY(N,SR,SI,XR,XI,INCX,YR,YI,INCY)                           
      DOUBLE PRECISION SR,SI,XR(1),XI(1),YR(1),YI(1),FLOP                       
      IF (N .LE. 0) RETURN                                                      
      IF (SR .EQ. 0.0D0 .AND. SI .EQ. 0.0D0) RETURN                             
      IX = 1                                                                    
      IY = 1                                                                    
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1                                       
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1                                       
      DO 10 I = 1, N                                                            
         YR(IY) = FLOP(YR(IY) + SR*XR(IX) - SI*XI(IX))                          
         YI(IY) = YI(IY) + SR*XI(IX) + SI*XR(IX)                                
         IF (YI(IY) .NE. 0.0D0) YI(IY) = FLOP(YI(IY))                           
         IX = IX + INCX                                                         
         IY = IY + INCY                                                         
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      DOUBLE PRECISION FUNCTION WDOTCR(N,XR,XI,INCX,YR,YI,INCY)                 
      DOUBLE PRECISION XR(1),XI(1),YR(1),YI(1),S,FLOP                           
      S = 0.0D0                                                                 
      IF (N .LE. 0) GO TO 20                                                    
      IX = 1                                                                    
      IY = 1                                                                    
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1                                       
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1                                       
      DO 10 I = 1, N                                                            
         S = FLOP(S + XR(IX)*YR(IY) + XI(IX)*YI(IY))                            
         IX = IX + INCX                                                         
         IY = IY + INCY                                                         
   10 CONTINUE                                                                  
   20 WDOTCR = S                                                                
      RETURN                                                                    
      END                                                                       
      DOUBLE PRECISION FUNCTION WDOTCI(N,XR,XI,INCX,YR,YI,INCY)                 
      DOUBLE PRECISION XR(1),XI(1),YR(1),YI(1),S,FLOP                           
      S = 0.0D0                                                                 
      IF (N .LE. 0) GO TO 20                                                    
      IX = 1                                                                    
      IY = 1                                                                    
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1                                       
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1                                       
      DO 10 I = 1, N                                                            
         S = S + XR(IX)*YI(IY) - XI(IX)*YR(IY)                                  
         IF (S .NE. 0.0D0) S = FLOP(S)                                          
         IX = IX + INCX                                                         
         IY = IY + INCY                                                         
   10 CONTINUE                                                                  
   20 WDOTCI = S                                                                
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WSCAL(N,SR,SI,XR,XI,INCX)                                      
      DOUBLE PRECISION SR,SI,XR(1),XI(1)                                        
      IF (N .LE. 0) RETURN                                                      
      IX = 1                                                                    
      DO 10 I = 1, N                                                            
         CALL WMUL(SR,SI,XR(IX),XI(IX),XR(IX),XI(IX))                           
         IX = IX + INCX                                                         
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WSWAP(N,XR,XI,INCX,YR,YI,INCY)                                 
      DOUBLE PRECISION XR(1),XI(1),YR(1),YI(1),T                                
      IF (N .LE. 0) RETURN                                                      
      IX = 1                                                                    
      IY = 1                                                                    
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1                                       
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1                                       
      DO 10 I = 1, N                                                            
         T = XR(IX)                                                             
         XR(IX) = YR(IY)                                                        
         YR(IY) = T                                                             
         T = XI(IX)                                                             
         XI(IX) = YI(IY)                                                        
         YI(IY) = T                                                             
         IX = IX + INCX                                                         
         IY = IY + INCY                                                         
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      DOUBLE PRECISION FUNCTION WNRM2(N,XR,XI,INCX)                             
      DOUBLE PRECISION XR(1),XI(1),PYTHAG,S                                     
C     NORM2(X)                                                                  
      S = 0.0D0                                                                 
      IF (N .LE. 0) GO TO 20                                                    
      IX = 1                                                                    
      DO 10 I = 1, N                                                            
         S = PYTHAG(S,XR(IX))                                                   
         S = PYTHAG(S,XI(IX))                                                   
         IX = IX + INCX                                                         
   10 CONTINUE                                                                  
   20 WNRM2 = S                                                                 
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RROTG(DA,DB,C,S)                                               
C                                                                               
C     CONSTRUCT GIVENS PLANE ROTATION.                                          
C                                                                               
      DOUBLE PRECISION DA,DB,C,S,RHO,PYTHAG,FLOP,R,Z                            
C                                                                               
      RHO = DB                                                                  
      IF ( DABS(DA) .GT. DABS(DB) ) RHO = DA                                    
      C = 1.0D0                                                                 
      S = 0.0D0                                                                 
      Z = 1.0D0                                                                 
      R = FLOP(DSIGN(PYTHAG(DA,DB),RHO))                                        
      IF (R .NE. 0.0D0) C = FLOP(DA/R)                                          
      IF (R .NE. 0.0D0) S = FLOP(DB/R)                                          
      IF ( DABS(DA) .GT. DABS(DB) ) Z = S                                       
      IF ( DABS(DB) .GE. DABS(DA) .AND. C .NE. 0.0D0 ) Z = FLOP(1.0D0/C)        
      DA = R                                                                    
      DB = Z                                                                    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WMUL(AR,AI,BR,BI,CR,CI)                                        
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI,T,FLOP                                 
C     C = A*B                                                                   
      T = AR*BI + AI*BR                                                         
      IF (T .NE. 0.0D0) T = FLOP(T)                                             
      CR = FLOP(AR*BR - AI*BI)                                                  
      CI = T                                                                    
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WDIV(AR,AI,BR,BI,CR,CI)                                        
      DOUBLE PRECISION AR,AI,BR,BI,CR,CI                                        
C     C = A/B                                                                   
      DOUBLE PRECISION S,D,ARS,AIS,BRS,BIS,FLOP                                 
      S = DABS(BR) + DABS(BI)                                                   
      IF (S .EQ. 0.0D0) CALL ZRROR(27)                                          
      IF (S .EQ. 0.0D0) RETURN                                                  
      ARS = AR/S                                                                
      AIS = AI/S                                                                
      BRS = BR/S                                                                
      BIS = BI/S                                                                
      D = BRS**2 + BIS**2                                                       
      CR = FLOP((ARS*BRS + AIS*BIS)/D)                                          
      CI = (AIS*BRS - ARS*BIS)/D                                                
      IF (CI .NE. 0.0D0) CI = FLOP(CI)                                          
      RETURN                                                                    
      END                                                                       
      SUBROUTINE WSIGN(XR,XI,YR,YI,ZR,ZI)                                       
      DOUBLE PRECISION XR,XI,YR,YI,ZR,ZI,PYTHAG,T                               
C     IF Y .NE. 0, Z = X*Y/ABS(Y)                                               
C     IF Y .EQ. 0, Z = X                                                        
      T = PYTHAG(YR,YI)                                                         
      ZR = XR                                                                   
      ZI = XI                                                                   
      IF (T .NE. 0.0D0) CALL WMUL(YR/T,YI/T,ZR,ZI,ZR,ZI)                        
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RROT(N,DX,INCX,DY,INCY,C,S)                                    
C                                                                               
C     APPLIES A PLANE ROTATION.                                                 
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S,FLOP                               
      INTEGER I,INCX,INCY,IX,IY,N                                               
C                                                                               
      IF (N.LE.0) RETURN                                                        
      IX = 1                                                                    
      IY = 1                                                                    
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1                                       
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1                                       
      DO 10 I = 1,N                                                             
        DTEMP = FLOP(C*DX(IX) + S*DY(IY))                                       
        DY(IY) = FLOP(C*DY(IY) - S*DX(IX))                                      
        DX(IX) = DTEMP                                                          
        IX = IX + INCX                                                          
        IY = IY + INCY                                                          
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE ZRROR(NUM)                                                     
      WRITE(6,6)NUM                                                             
    6 FORMAT(' SUBROUTINE ZRROR...',I5)                                         
      RETURN                                                                    
      END                                                                       
      DOUBLE PRECISION FUNCTION FLOP(X)                                         
      DOUBLE PRECISION X                                                        
C     SYSTEM DEPENDENT FUNCTION                                                 
C     COUNT AND POSSIBLY CHOP EACH FLOATING POINT OPERATION                     
C     FLP(1) IS FLOP COUNTER                                                    
C     FLP(2) IS NUMBER OF PLACES TO BE CHOPPED                                  
C                                                                               
      INTEGER SYM,SYN(4),BUF(256),CHAR,FLP(2),FIN,FUN,LHS,RHS,RAN(2)            
      COMMON /COM/ SYM,SYN,BUF,CHAR,FLP,FIN,FUN,LHS,RHS,RAN                     
C                                                                               
      DOUBLE PRECISION MASK(13),XX,MM                                           
      LOGICAL LX(2),LM(2)                                                       
      EQUIVALENCE (LX(1),XX),(LM(1),MM)                                         
      DATA MASK        / ZFFFFFFFFFFFFFFF0,ZFFFFFFFFFFFFFF00,                   
     $ ZFFFFFFFFFFFFF000,ZFFFFFFFFFFFF0000,ZFFFFFFFFFFF00000,                   
     $ ZFFFFFFFFFF000000,ZFFFFFFFFF0000000,ZFFFFFFFF00000000,                   
     $ ZFFFFFFF000000000,ZFFFFFF0000000000,ZFFFFF00000000000,                   
     $ ZFFFF000000000000,ZFFF0000000000000/                                     
C********                                                                       
      FLOP=X                                                                    
      RETURN                                                                    
C********                                                                       
C                                                                               
C     FLP(1) = FLP(1) + 1                                                       
C     K = FLP(2)                                                                
C     FLOP = X                                                                  
C     IF (K .LE. 0) RETURN                                                      
C     FLOP = 0.0D0                                                              
C     IF (K .GE. 14) RETURN                                                     
C     XX = X                                                                    
C     MM = MASK(K)                                                              
C     LX(1) = LX(1) .AND. LM(1)                                                 
C     LX(2) = LX(2) .AND. LM(2)                                                 
C     FLOP = XX                                                                 
C     RETURN                                                                    
      END                                                                       
      SUBROUTINE WRSCAL(N,S,XR,XI,INCX)                                         
      DOUBLE PRECISION S,XR(1),XI(1),FLOP                                       
      IF (N .LE. 0) RETURN                                                      
      IX = 1                                                                    
      DO 10 I = 1, N                                                            
         XR(IX) = FLOP(S*XR(IX))                                                
         IF (XI(IX) .NE. 0.0D0) XI(IX) = FLOP(S*XI(IX))                         
         IX = IX + INCX                                                         
   10 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      DOUBLE PRECISION FUNCTION PYTHAG(AR,AI)                                   
      DOUBLE PRECISION AR,AI                                                    
      PYTHAG=DSQRT(AR*AR+AI*AI)                                                 
      RETURN                                                                    
      END                                                                       














