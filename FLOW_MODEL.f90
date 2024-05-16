
    PROGRAM MAIN

    COMMON/FLOW1/NR,NC,NTS
    COMMON/FLOW2/DELX,DELY,DELT
    COMMON/FLOW3/S,W
    COMMON/FLOW4/TX,TY

    DIMENSION S(150,150),W(150,150,100),W_NMN(150,150,100),HIN(150,150),HK(150,150,100)
  	DIMENSION TX(150,150),TY(150,150),HF(150,150)
	DIMENSION IDB(1000),JDB(1000),INB(1000),JNB(1000)

	INTEGER NR,NC,NTS,NDRICHLETB
	REAL DELT,DELX,DELY,EPS

    OPEN (UNIT=11,FILE='INPUT.TXT',STATUS='OLD')
    READ (11,*) NR
    READ (11,*) NC
    READ (11,*) NTS
    READ (11,*) DELX
    READ (11,*) DELT
    READ (11,*) EPS

    CLOSE (11)

    DELY=DELX

    OPEN (UNIT=11,FILE='TX.TXT',STATUS='UNKNOWN')
	DO I=1,NR
    READ(11,*) (TX(I,J),J=1,NC-1)
    END DO
    CLOSE (11)

    OPEN (UNIT=11,FILE='TY.TXT',STATUS='UNKNOWN')
	DO J=1,NC
    READ(11,*) (TY(I,J),I=1,NR-1)
    END DO
    CLOSE (11)


    OPEN (UNIT=12,FILE='S.TXT',STATUS='UNKNOWN')
	DO I=1,NR
    READ(12,*) (S(I,J),J=1,NC)
    END DO
    CLOSE (12)


    OPEN (UNIT=13,FILE='W.TXT',STATUS='UNKNOWN')
    DO K=1,NTS
	DO I=1,NR
    READ(13,*) (W(I,J,K),J=1,NC)
    END DO
    END DO
    CLOSE (13)


    OPEN (UNIT=14,FILE='INITIAL_HEADS.TXT',STATUS='UNKNOWN')
	DO I=1,NR
    READ(14,*) (HIN(I,J),J=1,NC)
    END DO
    CLOSE (14)


    OPEN (UNIT=15,FILE='DIRICHILET_BOUNDARY_HEADS.TXT',STATUS='UNKNOWN')
    READ (15,*) N_DIRICHLET_BC
    DO I=1,N_DIRICHLET_BC
    READ(15,*) IDB(I),JDB(I),(HK(IDB(I),JDB(I),K),K=1,NTS)
    S(IDB(I),JDB(I))=1E34
    END DO
    CLOSE (15)


    OPEN (UNIT=16,FILE='NUEMAN_BOUNDARY_CONDITION.TXT',STATUS='UNKNOWN')
    READ (16,*) N_NEUMAN_BC
    DO I=1,N_NEUMAN_BC
    READ(16,*) INB(I),JNB(I),(W_NMN(INB(I),JNB(I),K),K=1,NTS)
    END DO
    CLOSE (16)


    DO K=1,NTS
	DO I=1,NR
    DO J=1,NC
    W(I,J,K)=W(I,J,K)+W_NMN(I,J,K)/DELX
    END DO
    END DO
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!C	SIMULATION STARTS
    OPEN (UNIT=14,FILE='FINAL_HEADS.TXT',STATUS='UNKNOWN')
	DO 100 K=1,NTS



    DO I=1,N_DIRICHLET_BC
	HIN(IDB(I),JDB(I))=HK(IDB(I),JDB(I),K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   WRITE (*,*) "K=",K
!	WRITE (*,*) IDB(I),JDB(I),HK(IDB(I),JDB(I),K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	END DO

	CALL SIMU(K,HIN,HF)

!	INITIALIZATION OF NEXT TIME STEP

		DO 1002 I=1,NR
		DO 1001	J=1,NC
        HIN(I,J)=HF(I,J)
1001	CONTINUE
1002	CONTINUE

        WRITE (14,*) "K=",K
        DO I=1,NR
        WRITE (14,14) (HF(I,J),J=1,NC)
        END DO

        14 FORMAT (10F8.2)
100     CONTINUE
    CLOSE (14)

    END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE SIMULATION

	SUBROUTINE SIMU(K,HIN,HF)

    COMMON/FLOW1/NR,NC,NTS
    COMMON/FLOW2/DELX,DELY,DELT
    COMMON/FLOW3/S,W
    COMMON/FLOW4/TX,TY

	DIMENSION S(150,150),W(150,150,100),HI(150,150),H(150,150),HIN(150,150)
   	DIMENSION A(150),B(150),C(150),D(150),HT(150),TX(150,150),TY(150,150)
	DIMENSION E(150),F(150),G(150),Y(150)
	DIMENSION DXX(150,150),DYY(150,150),HF(150,150),V(150,150)

	INTEGER NR,NC,K

	REAL DELX,DELY,EPS,DELT

	DO 22 I=1,NR
	DO 21 J=1,NC
    HI(I,J)=HIN(I,J)
21	CONTINUE
22	CONTINUE

!C	ROW WISE HEAD CALCULATION

1500	DO 90 I=1,NR
		DO 80 J=1,NC
			IF(J.EQ.1)THEN
				TX1=0.0
				DX1=0.00001
			ELSE
				TX1=TX(I,J-1)
				DX1=DELX
			ENDIF
			IF(I.EQ.1) THEN
				TY1=0.0
				DY1=0.00001
				H1=0.0
			ELSE
				TY1=TY(I-1,J)
				DY1=DELY
				H1=HI(I-1,J)
			ENDIF
			IF(J.EQ.NC) THEN
				TX2=0.0
				DX2=0.00001
			ELSE
				DX2=DELX
				TX2=TX(I,J)
			ENDIF
			IF(I.EQ.NR) THEN
				TY2=0.0
				DY2=0.00001
				H2=0.0

			ELSE
				TY2=TY(I,J)
				DY2=DELY
				H2=HI(I+1,J)
			ENDIF

			A(J)=(2*TX1)/((DX1+DX2)*DX1)
			B(J)= ((-2)*(TX2/DX2+TX1/DX1)/(DX1+DX2))-(2*(TY2/DY2+TY1/DY1)/(DY1+DY2))-(S(I,J)/DELT)
			C(J)= (2*TX2)/((DX1+DX2)*DX2)
			D(J)= ((-2)*((TY2*H2/DY2)+(TY1*H1/DY1))/(DY1+DY2))+W(I,J,K)-(S(I,J)*HIN(I,J)/DELT)
80		CONTINUE
			 N=NC
     		CALL THOMALGOL(N,A,B,C,D,HT)

		DO 91 J=1,NC
			H(I,J)=HT(J)
91		CONTINUE

90	CONTINUE

!C	COLUMN WISE CALCULATION

	DO 182 J=1,NC
		DO 181 I=1,NR
			IF(J.EQ.1) THEN
				TX1=0.0
				DX1=0.00001
				H4=0.0
			ELSE
				TX1=TX(I,J-1)
				DX1=DELX
				H4=H(I,J-1)
			ENDIF
			IF(I.EQ.1) THEN
				TY1=0.0
				DY1=0.00001
			ELSE
				TY1=TY(I-1,J)
				DY1=DELY
			ENDIF
			IF(J.EQ.NC) THEN
				TX2=0.0
				DX2=0.00001
				H3=0.0
			ELSE
				DX2=DELX
				TX2=TX(I,J)
				H3=H(I,J+1)
			ENDIF
			IF(I.EQ.NR) THEN
				TY2=0.0
				DY2=0.00001
			ELSE
				TY2=TY(I,J)
				DY2=DELY
			ENDIF
			A(I)= 2*TY1/((DY1+DY2)*DY1)
			B(I)= ((-2)*(TX2/DX2+TX1/DX1)/(DX1+DX2))-(2*(TY2/DY2+TY1/DY1)/(DY1+DY2))-(S(I,J)/DELT)
			C(I)= 2*TY2/((DY1+DY2)*DY2)
			D(I)=((-2)*((TX2*H3/DX2)+(TX1*H4/DX1))/(DX1+DX2))+W(I,J,K)-((S(I,J)*HIN(I,J))/DELT)

181		CONTINUE
		N=NR

	    CALL THOMALGOL(N,A,B,C,D,HT)

		DO 192 I=1,NR
			HF(I,J)=HT(I)
192		CONTINUE


182	CONTINUE


!C	CONVERGENCE CRITERA
        DIFF=0.0
		DO 220 I=1,NR
        DO 210 J=1,NC

        TEMP= ABS(HF(I,J)-H(I,J))

        IF(TEMP.GT.DIFF) THEN
        DIFF=TEMP
        ENDIF

210		CONTINUE
220		CONTINUE
		IF(DIFF.LT.EPS) THEN
        DO 240 I=1,NR
        DO 230 J=1,NC
            HI(I,J)=HF(I,J)
230		CONTINUE
240		CONTINUE
        ITERATION=ITERATION+1
        GO TO 1500
		ELSE

		ENDIF

	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE THOMALGOL   (THOMAS ALGORITHM)

	SUBROUTINE THOMALGOL(N,A,B,C,D,HT)

	DIMENSION CSTAR(150),DSTAR(150),HT(150)
	DIMENSION A(150),B(150),C(150),D(150)

	DO 53 J=1,N-1
		IF(J.EQ.1) THEN
			CSTAR(J)=C(J)/B(J)
		ELSE
			CSTAR(J)=C(J)/(B(J)-(CSTAR(J-1)*A(J)))
		ENDIF
53	CONTINUE

	DO 55 J=1,N
		IF(J.EQ.1) THEN
			DSTAR(J)=D(J)/B(J)
		ELSE
			DSTAR(J)=(D(J)-(DSTAR(J-1)*A(J)))/(B(J)-(CSTAR(J-1)*A(J)))
		ENDIF
55	CONTINUE
	HT(N)=DSTAR(N)

	DO 502 J=2,N
		JJ=N-J+1
		HT(JJ)= DSTAR(JJ)-(HT(JJ+1)*CSTAR(JJ))
502	CONTINUE
	RETURN
	END



