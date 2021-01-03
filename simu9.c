/*
 * initialize  [   0]
 * hot to cold [   1]-[ 500]
 * mag         [ 501]-[1000]
 * cold to hot [1001]-[1500]
 * demag       [1501]-[2000]
 *
 * ������ñ��bed
 * h��m_dot��¸
 * ���ʺ�ʬˡ

 * AMR�٥å�¿�ز����륽����(2008/4/21����)
���ĥ������ȹ礦������ǧ�Ѥߡ�
�߼�Ψ������(M_PI or 3.14)

�ⲹü����layer1,2,,,���������롣

��ư��He��ʪ���ͤ��׻�����simu9
simu8����vis�η夬1��¿����
ref,COP�ʤɸ������Ф�0.01%�ʲ�

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define SPACE           (200)    //space step 200
#define TIME           (2001)    //time  step 2001
#define N_2            (1000) //1000
#define N_4             (500) //500

#define NUM            (1500) //1000
#define CYCLE_END      (1500) //1000
#define TAU             (8.0)    //1cycle [sec] 8.0
#define END          (1.0e-4)    //1.0e-4

	double Treg[SPACE][TIME];	//磁性体の温度を格納する二次元配列
	double Tgas[SPACE][TIME]; //流体の温度を格納する二次元配列

/************************ AUTO m_dot mode *******************/
#define AUTO
#ifdef AUTO
#define M_START     (2.0e-3)	//熱交換流体の質量の下限値
#define M_END       (3.0e-2)	//熱交換流体の質量の上限値
#define DM          (5.0e-4)	//熱交換流体の質量の変化幅
#endif

#define T_HIGH         (36.0)    //high end[K]
#define T_LOW          (31.0)    //low  end[K]
#define GAS_FLOW     (1.0e-3)    //m_dot[kg/sec](1秒間に流れる熱交換流体の質量)

//-- LAYER parameter -------------------------------------
#define LAYER             (2)    //層数（例：2層の場合→2) ※

int layer_parameter[LAYER+1];	//層の区切り目+両端

//-- magnet density [kg/m3] --------------------------------
void init_layer(double rho_r[LAYER], int layer_parameter[LAYER+1]);
/*********************/
/* RAl2    (6.10e3)  */
/* ErCo2   (1.03e4)  */
/* GSG    (7.955e3)  */
/*********************/
//#define RHO_MAG1      (6.10e3)    //RAl2
//#define RHO_MAG1      (1.03e4)    //ErCo2
//#define RHO_MAG1     (7.955e3)    //GSG

/* 各磁性体の密度、 充填する位置を与える関数 ※ */
void init_layer(double rho_r[LAYER], int layer_parameter[LAYER+1])
{
	rho_r[0] = 10.3e3;	//1層目の磁性体の密度
	rho_r[1] = 9.73e3;	//2層目の磁性体の密度
	layer_parameter[0] = 0;		//高温端の位置
	layer_parameter[1] = 200;	//層の区切りの位置
	layer_parameter[2] = 200;	//低温端の位置
}


//== prototype declaration ===============================================
void initialize(double Tgas[SPACE][TIME], double Treg[SPACE][TIME]);
void hot_to_cold(double Thigh, double fac1,
									double fac2[LAYER], double data_C[LAYER][3][NUM],
									double Tgas[SPACE][TIME], double Treg[SPACE][TIME]);
void cold_to_hot(double Tlow, double fac1, double fac2[LAYER],
									double data_C[LAYER][3][NUM],
									double Tgas[SPACE][TIME], double Treg[SPACE][TIME]);
void mag_demag(int time, double Treg[SPACE][TIME], double data_S[LAYER][3][NUM]);
double capacity(int lay, int field, double data_C[LAYER][3][NUM], double imput_temp);
double cal_angle(double data1_A, double data1_B, double data1_C,
								 double data2_A, double data2_B);
int end_check(int cycle, double Tgas[SPACE][TIME]);
void output(FILE *fp_Tgas, FILE *fp_Treg, FILE *fp_simu9_auto, FILE *fp_simu9_data,
						int cycle, double m_dot, double tau_flow,
						double dt, double C_gas, double Ntu, double lambda[LAYER],
						double Treg[SPACE][TIME], double Tgas[SPACE][TIME],
						double h, double rho_r[LAYER], double gas_data[5]);
void gnuplot_out(double Treg[SPACE][TIME]);
int layer_change(int step, int layer_parameter[LAYER+1]);

double capacity_old(int column, double data_C[5][15], double imput_temp);

//--------------------------------------------------------------------------//

int main(void)
{
	int i, j;
	int cycle;
	int check;
	int n2, n4;
	int n_capa, n_entropy, n_He;
	int lay;
	double fac1, fac2[LAYER];
	double L,r,Ac,dh,eps,as,h,Aw;
	double rho_r[LAYER];
	double dx, dt;
	double Thigh, Tlow, Tmid;
	double tau_cycle, tau_flow;
	double Ntu, lambda[LAYER];
	double m_dot;
	double C_gas;
	double rho_f;
	double vis_gas;
	double k_gas;
	double M_reg[LAYER];
	double Re, Pr, Nu;  //non dimension NUM
	double data_C[LAYER][3][NUM]; /*[material1,2][Temp,0T,5T][NUM]*/
	double data_S[LAYER][3][NUM]; /*[material1,2][Temp,0T,5T][NUM]*/
//	int layer_parameter[LAYER+1];
	double He5atm[5][15];
//Temp[K]	density[kg/m3]	Cp[J/kg-K]	viscosity[kg/m-s]	k[W/m-K]=[J/m-K-s]
	double gas_data[5];

	//== file open ==================================================
	FILE *fp_Tgas, *fp_Treg;
	FILE *fp_capa[LAYER];
	FILE *fp_entropy[LAYER];
	FILE *fp_simu9_auto;
	FILE *fp_simu9_data;
	FILE *fp_simu9_He5atm;

	/* シミュレーション結果書き込み用ファイルの作成 */
	fp_Tgas         = fopen("Tgas9.txt",       "w");
	fp_Treg         = fopen("Treg9.txt",       "w");
	fp_simu9_auto   = fopen("simu9-out.txt",   "w");
	fp_simu9_data   = fopen("simu9-data.txt",  "w");

	/* ファイルの読み込み */
	fp_simu9_He5atm = fopen("He5atm-data.txt", "r");
	fp_capa[0]      = fopen("(ErDy)Co2_kg/ErCo2_C_kg.txt", "r"); /* ※ */
	fp_capa[1]      = fopen("LaNi2_kg_GdNi2_select/C_Tc=70K_Gd3+_LaNi2.txt", "r");	/* ※ */
	fp_entropy[0] = fopen("(ErDy)Co2_kg/ErCo2_S_kg.txt", "r");	/* ※ */
	fp_entropy[1] = fopen("LaNi2_kg_GdNi2_select/S_Tc=70K_Gd3+_LaNi2.txt", "r");	/* ※ */

/*
C_ErCo2_0T_5T_02.txt
S_ErCo2_0T_5T.txt
C_HoAl2_0T_5T_hokan.txt
S_HoAl2_0T_5T_hokan.txt
C_DyAl2_0T_5T_trace.txt
S_DyAl2_0T_5T_trace.txt
C_kasou_GdSiGe_base.txt
S_kasou_GdSiGe_base.txt
material9/C_DyAl2_0T_5T_trace.txt
material9/S_DyAl2_0T_5T_trace.txt
material9/C_HoAl2_0T_5T_hokan.txt
material9/S_HoAl2_0T_5T_hokan.txt
debyesample9/C_DyAl2_20K_debye.txt
debyesample9/S_DyAl2_20K_debye.txt
*/

/*ファイルが用意されているかどうかの確認 */
	if(fp_Tgas == NULL){
		printf("can't open fpTgas !\n");
		exit(1);
	}
	if(fp_Treg == NULL){
		printf("can't open fpTreg !\n");
		exit(1);
	}
	if(fp_simu9_auto == NULL){
		printf("can't open fp_simu9_auto !!\n");
		exit(1);
	}
	if(fp_simu9_data == NULL){
		printf("can't open fp_simu9_data !!\n");
		exit(1);
	}
	for(i = 0; i < LAYER; i++){
		if(fp_capa[i] == NULL){
			printf("can't open fp_capa[%d] !!\n", i);
			exit(1);
		}
		if(fp_entropy[i] == NULL){
			printf("can't open fp_entropy[%d] !!\n", i);
			exit(1);
		}
	}
	if(fp_simu9_He5atm == NULL){
		printf("can' open fp_simu9_He5atm !!\n");
		exit(1);
	}

	/* ファイルのヘッダー作成 */
	fprintf(fp_simu9_auto, "Thigh Tlow m_dot(kg/sec) ref_load(J) warm_load COP FOM h cycle Ntu\n");
	fprintf(fp_Tgas, "TimeStep coldend hotend m_ dot\n");
	fprintf(fp_Treg, "SpaceStep hot-to-cold mag cold-to-hot demag m_dot cycle\n");


	for(i = 0; i < LAYER; i++){
		n_capa = 0;
		while(EOF != (fscanf(fp_capa[i], "%lf %lf %lf",
						&data_C[i][0][n_capa],&data_C[i][1][n_capa],&data_C[i][2][n_capa]))){
			n_capa++;
		}

		n_entropy = 0;
		while(EOF != (fscanf(fp_entropy[i], "%lf %lf %lf",
					&data_S[i][0][n_entropy],&data_S[i][1][n_entropy],&data_S[i][2][n_entropy]))){
			n_entropy++;
		}
	}

	n_He = 0;
	while(EOF != (fscanf(fp_simu9_He5atm, "%lf %lf %lf %lf %lf",
				&He5atm[0][n_He],&He5atm[1][n_He],&He5atm[2][n_He],&He5atm[3][n_He],&He5atm[4][n_He]))){
		n_He++;
	}


	//-----------------------------------------------------------//
	n2        = N_2;//(int)(TIME/2);
	n4        = N_4;//(int)(TIME/4);
	tau_cycle = TAU;            		//[sec]:1cycle
	tau_flow  = tau_cycle/4.0;  		//[sec]:1flow
	Thigh     = T_HIGH;         		//[K]
	Tlow      = T_LOW;          		//[K]
	Tmid      = 0.50 * (Thigh + Tlow);	//(K)
	m_dot     = GAS_FLOW;       		//[kg/sec]

	//== bed parameters ==============
	L         = 0.1;            //[m]:AMRベッドの長さ
	r         = 0.0225;         //[m]:AMRベッド底面の半径
	Ac        = M_PI*r*r;       //[m2]:AMRベッド底面の断面積
//	Ac        = 1.59e-3; /*ikeda*/

	//== magnet parameters ===========
	dh        = 4.0e-4;         //[m]:磁性体粒径（直径）
//	rho_r[0]  = RHO_MAG1;        //[kg/m3]:magnetic density
	init_layer(rho_r, layer_parameter);
	eps       = 0.36;           //AMRベッドの空孔率
	as        = 6.0*(1.0 - eps)/dh; //単位体積あたりの磁性体-流体間の熱交換面積（比表面積 [m^2/m^3]）

	//== gas parameters (Helium) =====
	for(i = 1; i < 5; i++){
		gas_data[i] = capacity_old(i, He5atm, Tmid);
	}
	rho_f     = gas_data[1];	//[kg/m3] :流体(Heガス)の密度
	C_gas     = gas_data[2];	//[J/kg-K]:流体の比熱
	vis_gas   = gas_data[3];	//[kg/m-sec]:流体の粘性係数
	k_gas     = gas_data[4];	//[J/m-K-sec]:流体の熱伝導係数
	h         = 0.0;      		//[W/m2-K] = [J/sec-m2-K]:heat conduction ratio

	//== other parameters ============
	Aw        = Ac*as;          //[]

	for(lay = 0; lay < LAYER; lay++){
		M_reg[lay] = Ac*L*(1.0 - eps)*rho_r[lay]; //[kg]:AMRベッド内に充填している各磁性体の質量
	}


	dx = L/(1.0*SPACE - 1.0);	//[m]:空間の分割幅
	dt = tau_cycle/(1.0*TIME - 1.0); //[sec]:時間の分割幅

	/* シミュレーション条件を出力 */
	fprintf(stderr, "高温端温度 = %g [K]\n", Thigh);
	fprintf(stderr, "低温端温度 = %g [K]\n", Tlow);
	fprintf(stderr, "中間温度 = %g [K]\n", Tmid);
	fprintf(stderr, "流体の密度 = %g [kg/m^3]\n", rho_f);
	fprintf(stderr, "流体の比熱 = %g [J/kg-K]\n", C_gas);
	fprintf(stderr, "流体の粘性係数 = %g [kg/m-s]\n", vis_gas);
	fprintf(stderr, "流体の熱伝導係数  = %g [J/m-K-s]\n", k_gas);

	for(lay = 0; lay < LAYER; lay++){
		fprintf(stderr, "%d層目磁性体密度 = %g [kg/m^3]\n",lay+1, rho_r[lay]);
		fprintf(stderr, "%d層目磁性体質量 = %g [kg]\n",lay+1, M_reg[lay]);
	}
	for(lay = 0;lay <= LAYER; lay++){
		fprintf(stderr, "layer_parameter = %d\n",layer_parameter[lay]);
	}
	fprintf(stderr, "dx = %g [m]\n", dx);
	fprintf(stderr, "dt = %g [s]\n", dt);
	fprintf(stderr, "比表面積 = %g [1/m]\n", as);
	fprintf(stderr, "AMRベッド断面積 = %g [m^2]\n", Ac);
	fprintf(stderr, "Aw = %g [m]\n", Aw);
	fprintf(stderr, "m_dot  = %g [kg/s]\n", m_dot);
	fprintf(stderr, "tau_cycle = %g [s]\n", tau_cycle);
	fprintf(stderr, "tau_flow  = %g [s]\n\n", tau_flow);

	#ifdef AUTO
	fprintf(stderr, "RUNNING AUTO MODE !!\n");
	clock_t start_all, end_all;
	start_all = clock();
	/* M_STARTの値からDMずつ足していき、M_ENDよりm_dotが大きくなるまで処理を繰り返す */
	for(m_dot = M_START; m_dot <= M_END; m_dot += DM){
		clock_t start, end;
		start = clock();
	#endif

		Re = (m_dot*dh)/(Ac*vis_gas); //レイノルズ数の導出
		Pr = (vis_gas*C_gas)/k_gas;
		Nu = 2.0 + 1.1*pow(Re, 0.6)*pow(Pr, 0.33333);
		h  = (Nu*k_gas)/dh;
//		h = 3000.0;

		fprintf(stderr, "m_dot = %g h = %f\n", m_dot, h);

		Ntu = (h*Aw*L)/(m_dot*C_gas);
		fac1 = (dx * h * Aw)/(m_dot * C_gas);

		for(lay = 0; lay < LAYER; lay++){
			/*--1.0/Gamma = lambda = (m_dot*C_gas*tau_flow)/(Mreg*C_reg)--*/
			lambda[lay] = (m_dot*C_gas*tau_flow)/(M_reg[lay]); //omit (1.0/C_reg)

			fac2[lay] = (dt * h * Aw * L)/(M_reg[lay]); //1.0/C_reg
		}

		initialize(Tgas, Treg);

		//== main loop ==============================================
		for(cycle = 0;  cycle < CYCLE_END; cycle++){
//			printf("cycle = %d\n", cycle);
			hot_to_cold(Thigh, fac1, fac2, data_C, Tgas, Treg);
			mag_demag(n4+1, Treg, data_S);
			cold_to_hot(Tlow, fac1, fac2, data_C, Tgas, Treg);
			mag_demag(n4+n2+1, Treg, data_S);
			check = end_check(cycle, Tgas);
			if(check) break; //convergence OK:check == 1, NO:check == 0
		}

		output(fp_Tgas, fp_Treg, fp_simu9_auto, fp_simu9_data,
						cycle, m_dot, tau_flow,
						dt, C_gas, Ntu, lambda,
						Treg, Tgas, h, rho_r, gas_data);

	#ifdef AUTO
		end = clock();
		printf("time = %.2f秒\n", (double)(end-start)/CLOCKS_PER_SEC);
	}
	#endif

	//gnuplot_out(Tgas);

	/*
	for(i = 0; i < SPACE; i++){
		printf("%g %g %g %g %g %g %g %g \n",
						Treg[i][n2+n4], Treg[i][n2+n4+1], Treg[i][TIME-1], Treg[i][0],
						Tgas[i][n2+n4], Tgas[i][n2+n4+1], Tgas[i][TIME-1], Tgas[i][0]);
	}
	*/

	fclose(fp_Tgas);
	fclose(fp_Treg);
	fclose(fp_capa[LAYER-1]);
	fclose(fp_capa[LAYER-2]);
	fclose(fp_entropy[LAYER-1]);
	fclose(fp_entropy[LAYER-2]);
	fclose(fp_simu9_auto);
	fclose(fp_simu9_data);
	fclose(fp_simu9_He5atm);

	end_all = clock();
	printf("プログラム全体の処理時間 = %.2f秒\n", (double)(end_all-start_all)/CLOCKS_PER_SEC);
	return 0;
}

void initialize(double Tgas[SPACE][TIME], double Treg[SPACE][TIME])
{
	int i, j;
	double deltaTemp;
	double Thigh = T_HIGH;
	double Tlow  = T_LOW;

	for(j = 0; j < TIME; j++){
		for(i = 0; i < SPACE; i++){
			Tgas[i][j] = 0.0;
			Treg[i][j] = 0.0;
		}
	}
	deltaTemp = (Thigh - Tlow)/(1.0 * SPACE - 1.0);
	for(i = 0; i < SPACE; i++){
		Tgas[i][0] = Thigh - 1.0*i*deltaTemp;
		Treg[i][0] = Thigh - 1.0*i*deltaTemp;
	}
}

void hot_to_cold(double Thigh, double fac1,
									double fac2[LAYER], double data_C[LAYER][3][NUM],
									double Tgas[SPACE][TIME], double Treg[SPACE][TIME])
{
	int i, j;
	int lay;
	int field0T = 1;
	//int n2 = N_2;
	int n4 = N_4;
	double C_reg;

	for(j = 1; j <= n4; j++){
//		printf("j = %d\n",j);
		Tgas[0][j] = Thigh;
		for(i = 0; i < SPACE; i++){
			lay = layer_change(i,layer_parameter);
			C_reg = capacity(lay, field0T, data_C, Treg[i][j-1]);
			Treg[i][j] = Treg[i][j-1] + (fac2[lay]/C_reg)*(Tgas[i][j-1] - Treg[i][j-1]);
		}
		for(i = 0; i < SPACE-1; i++){
			Tgas[i+1][j] = Tgas[i][j] + fac1*(Treg[i][j] - Tgas[i][j]);
		}
	}
}


void cold_to_hot(double Tlow, double fac1,
									double fac2[LAYER], double data_C[LAYER][3][NUM],
									double Tgas[SPACE][TIME], double Treg[SPACE][TIME])
{
	int i, j;
	int lay;
	int field5T = 2;
	int n2 = N_2;
	int n4 = N_4;
	double C_reg;

	for(j = n2+1; j <= n2+n4; j++){
		Tgas[SPACE-1][j] = Tlow;
		for(i = SPACE-1; i >= 0; i--){
			lay = layer_change(i,layer_parameter);
			C_reg = capacity(lay, field5T, data_C, Treg[i][j-1]);
			Treg[i][j] = Treg[i][j-1] + (fac2[lay]/C_reg)*(Tgas[i][j-1] - Treg[i][j-1]);
		}
		for(i = SPACE-1; i > 0; i--){
			Tgas[i-1][j] = Tgas[i][j] + fac1*(Treg[i][j] - Tgas[i][j]);
		}
	}
}


/********************************************/
/*  imput Temp_C(0T) --> return Temp_C(5T)  */
/*                                          */
/*                                          */
/*          |  T   | S(0T)     | S(5T)    | */
/*          |------|-----------|----------| */
/*          |Temp_A|Entrpoy_A  |          | */
/*Temp_C -> |     --> Entropy_C|          | */
/*          |Temp_B|Entropy_B  |          | */
/*          |      |           |          | */
/*      <---|Temp_C| <-----------Entropy_C| */
/********************************************/

void mag_demag(int time, double Treg[SPACE][TIME], double data_S[LAYER][3][NUM])
{
	int i, j;
	int n2 = N_2;
	int n4 = N_4;//(int)(TIME/4);

	int before, after; // 0T to 5T or 5T to 0T
	int lay;

	double entropy_A;
	double entropy_B;
	double entropy_C;
	double temp_A;
	double temp_B;
	double temp_C;/*imput temp*/

	if(time == n4+1){
		before = 1; //0T
		after  = 2; //5T
	}else if(time == n4+n2+1){
		before = 2; //5T
		after  = 1; //0T
	}else{
		printf("err mag demag!!\n");
		exit(1);
	}

	for(i = 0; i < SPACE; i++){
		temp_C = Treg[i][time-1];

		lay = layer_change(i,layer_parameter);

		//printf("temp_C=%g\n", temp_C);
		for(j = 0; j < NUM; j++){
			if(temp_C < data_S[lay][0][j]){
				temp_A    = data_S[lay][0][j-1];
				temp_B    = data_S[lay][0][j];
				entropy_A = data_S[lay][before][j-1];
				entropy_B = data_S[lay][before][j];
				break;
			}
		}
		if((j == NUM)||(j == 0)){
			printf("err mag_demag1!!\n");
			exit(1);
		}
		entropy_C = cal_angle(temp_A, temp_B, temp_C, entropy_A, entropy_B);

		for(j = 0; j < NUM; j++){
			if(entropy_C < data_S[lay][after][j]){
				temp_A    = data_S[lay][0][j-1];
				temp_B    = data_S[lay][0][j];
				entropy_A = data_S[lay][after][j-1];
				entropy_B = data_S[lay][after][j];
				break;
			}
		}
		if((j == NUM)||(j == 0)){
			printf("err mag_demag2!!\n");
			exit(1);
		}
		temp_C = cal_angle(entropy_A, entropy_B, entropy_C, temp_A, temp_B);

		for(j = 0; j < n4; j++){
			Treg[i][time+j] = temp_C;
			Tgas[i][time+j] = temp_C;
		}
		if(time == n2+n4+1){
			Treg[i][0] = temp_C;
			Tgas[i][0] = temp_C;
		}
	}
}


double cal_angle(double data1_A, double data1_B, double data1_C,
								 double data2_A, double data2_B)
{
	double angle;
	double delta;
	double data2_C;

	angle   = (data1_C - data1_A)/(data1_B - data1_A);
	delta   = data2_B - data2_A;
	data2_C = angle * delta + data2_A;
	return data2_C;
}


/********************************************/
/*  imput Temp_C(0T) --> return Capa_C(0T)  */
/*                                          */
/*                                          */
/*          |  T   | C(0T)    | C(5T) |     */
/*          |------|----------|-------|     */
/*          |Temp_A|          |       |     */
/*Temp_C -> |     -|-> Capa_C-|-------|->   */
/*          |Temp_B|          |       |     */
/*          |      |          |       |     */
/*          |      |          |       |     */
/********************************************/

double capacity(int lay, int field, double data_C[LAYER][3][NUM], double imput_temp)
{
	int i;
	double capa_A;
	double capa_B;
	double capa_C;
	double temp_A;
	double temp_B;
	double temp_C = imput_temp;

	for(i = 0; i < NUM; i++){
		if(temp_C < data_C[lay][0][i]){
			 temp_A = data_C[lay][0][i-1];
			 temp_B = data_C[lay][0][i];
			 capa_A = data_C[lay][field][i-1];
			 capa_B = data_C[lay][field][i];
			break;
		}

		if(i == NUM-1){
			printf("err capacity !!\n");
			printf("tempC = %f, i = %d\n", temp_C, i);
			exit(1);
		}
	}
	capa_C = cal_angle(temp_A, temp_B, temp_C, capa_A, capa_B);
	return capa_C;
}


int end_check(int cycle, double Tgas[SPACE][TIME])
{
	int n2 = N_2;
	int n4 = N_4;
	double Tgas_end[2][CYCLE_END];

	Tgas_end[0][cycle] = Tgas[SPACE-1][n4]; //[0]is low end(hot to cold)
	Tgas_end[1][cycle] = Tgas[0][n2+n4];    //[1]is high end(cold to hot)

	if(cycle > 1){
		if((fabs(Tgas_end[0][cycle] - Tgas_end[0][cycle-1]) < END) &&
			 (fabs(Tgas_end[1][cycle] - Tgas_end[1][cycle-1]) < END)){
			return 1;
		}
	}
	return 0;
}

void output(FILE *fp_Tgas, FILE *fp_Treg, FILE *fp_simu9_auto, FILE *fp_simu9_data,
						int cycle, double m_dot, double tau_flow,
						double dt, double C_gas, double Ntu, double lambda[LAYER],
						double Treg[SPACE][TIME], double Tgas[SPACE][TIME],
						double h, double rho_r[LAYER], double gas_data[5])
{
	int i, j;
	int n2 = N_2;
	int n4 = N_4;
	int lay;
	double Thigh = T_HIGH;
	double Tlow  = T_LOW;
	double H_high, H_low;
	double ref_load, warm_load;
	double Tgas_high_end_sum = 0.0;
	double Tgas_low_end_sum  = 0.0;
	double COP, COPcar, FOM;
/*
	double Tgas_high_end_eff;
	double Tgas_low_end_eff;
	double thermal_eff1;
	double thermal_eff2;
	double lambda_0T_high, lambda_0T_low;
	double lambda_5T_high, lambda_5T_low;
*/
	double C_reg_0T_high, C_reg_0T_low; //0T_capa high & low end
	double C_reg_5T_high, C_reg_5T_low; //5T_capa high & low end

/*
	C_reg_0T_high = capacity(data_C_Temp, data_C_0T, Thigh);
	C_reg_0T_low  = capacity(data_C_Temp, data_C_0T, Tlow );
	C_reg_5T_high = capacity(data_C_Temp, data_C_5T, Thigh);
	C_reg_5T_low  = capacity(data_C_Temp, data_C_5T, Tlow );

	lambda_0T_high = lambda/C_reg_0T_high;
	lambda_0T_low  = lambda/C_reg_0T_low;
	lambda_5T_high = lambda/C_reg_5T_high;
	lambda_5T_low  = lambda/C_reg_5T_low;
*/
	for(j = 1; j <= n4; j++){
		Tgas_low_end_sum  += Tgas[SPACE-1][j];
		Tgas_high_end_sum += Tgas[0][n2+j];
	}
/*
	//output gas average tempetature
	Tgas_low_end_eff  = Tgas_low_end_sum /((double)(n4));
	Tgas_high_end_eff = Tgas_high_end_sum/((double)(n4));
	//printf("%g %g\n", Tgas_low_end_eff, Tgas_high_end_eff);

	//non dimension parameters
	thermal_eff1 = (Thigh - Tgas_low_end_eff)/(Thigh - Tlow);
	thermal_eff2 = (Tgas_high_end_eff - Tlow)/(Thigh - Tlow);
*/
	//discharge enthalpy
	H_low  = m_dot * dt * C_gas * Tgas_low_end_sum;
	H_high = m_dot * dt * C_gas * Tgas_high_end_sum;

	//refrigeration & warming
	ref_load  = (m_dot * Tlow * C_gas * tau_flow) - H_low;
	warm_load = H_high - (m_dot * Thigh * C_gas * tau_flow);

	//Cofficient of Performance
	COP    = ref_load/(warm_load - ref_load);
	COPcar = 1.0*T_LOW/(T_HIGH - T_LOW);
	FOM    = COP/COPcar;


	#ifndef AUTO
	fprintf(fp_simu9_data, "===== File : %s ==== Date : %s %s =====\n",
												 __FILE__, __DATE__, __TIME__);
	fprintf(fp_simu9_data, "Thigh     = %g(K)\n", Thigh);
	fprintf(fp_simu9_data, "Tlow      = %g(K) [Tmid = %g(K)]\n", Tlow, 0.5*(Thigh+Tlow));
	fprintf(fp_simu9_data, "m_dot     = %f(kg/sec)\n", m_dot);
	fprintf(fp_simu9_data, "ref_load  = %f(J)\n", ref_load);
	fprintf(fp_simu9_data, "warm_load = %f(J)\n", warm_load);
	fprintf(fp_simu9_data, "COP    = %f\n", COP);
	fprintf(fp_simu9_data, "FOM    = %f\n", FOM);
	fprintf(fp_simu9_data, "cycle  = %d\n", cycle);
	fprintf(fp_simu9_data, "Ntu    = %f\n", Ntu);
/*
	fprintf(fp_simu9_data, "lambda = %f ([0T]Thigh)\n", lambda_0T_high);
	fprintf(fp_simu9_data, "lambda = %f ([0T]Tlow) \n", lambda_0T_low );
	fprintf(fp_simu9_data, "lambda = %f ([5T]Thigh)\n", lambda_5T_high);
	fprintf(fp_simu9_data, "lambda = %f ([5T]Tlow) \n", lambda_5T_low );
	fprintf(fp_simu9_data, "thermal_eff1 = %f (Thigh-eff)\n", thermal_eff1);
	fprintf(fp_simu9_data, "thermal_eff2 = %f (eff-Tlow)\n",  thermal_eff2);

	fprintf(fp_simu9_data, "C_gas   = %e(J/kg-K)\n"    , C_Helium);
	fprintf(fp_simu9_data, "rho_gas = %e(kg/m3)\n"     , Density_Helium);
	fprintf(fp_simu9_data, "vis_gas = %e(kg/m-sec)\n"  , Viscosity_Helium);
	fprintf(fp_simu9_data, "k_gas   = %e(J/sec-K-m)\n" , k_Helium);
*/
/*
	rho_f     = gas_data[1];	//[kg/m3] :fluid density
	C_gas     = gas_data[2];	//[J/kg-K]:fluid specific heat
	vis_gas   = gas_data[3];	//[kg/m-sec]
	k_gas     = gas_data[4];	//[J/m-K-sec]
*/
	fprintf(fp_simu9_data, "C_gas   = %e(J/kg-K)\n"    , gas_data[2]);
	fprintf(fp_simu9_data, "rho_gas = %e(kg/m3)\n"     , gas_data[1]);
	fprintf(fp_simu9_data, "vis_gas = %e(kg/m-sec)\n"  , gas_data[3]);
	fprintf(fp_simu9_data, "k_gas   = %e(J/sec-K-m)\n" , gas_data[4]);

	fprintf(fp_simu9_data, "h       = %e(J/sec-K-m2)\n", h);
	for(lay = 0; lay < LAYER; lay++){
		fprintf(fp_simu9_data, "rho_r[%d]  = %f\n", lay, rho_r[lay]);
	}
	for(lay = 0;lay <= LAYER; lay++){
		fprintf(fp_simu9_data, "layer_parameter = %d\n",layer_parameter[lay]);
	}
	#endif

/*
	fprintf(fp_simu9_auto, "%f %f %f %f %f %f %f %f %d %f %f %f %f %f\n",
													Thigh, Tlow,
													m_dot, ref_load, warm_load, COP, FOM, h,
													cycle, Ntu,
													lambda_0T_high, lambda_0T_low,
													lambda_5T_high, lambda_5T_low);
*/
	fprintf(fp_simu9_auto, "%f %f %f %f %f %f %f %f %d %f\n",
													Thigh, Tlow,
													m_dot, ref_load, warm_load, COP, FOM, h,
													cycle, Ntu);


	// #ifndef AUTO
	//high_end & low_end output gas
	for(j = 1; j <= n4; j++){

		fprintf(fp_Tgas, "%d %f %f\n", j, Tgas[SPACE-1][j], Tgas[0][n2+j]);
	}
	//High to Low & L to H last time temperature
	for(i = 0; i < SPACE; i++){
		fprintf(fp_Treg, "%d %f %f %f %f %f\n", i,
											Treg[i][n4], Treg[i][n2],
											Treg[i][n2+n4], Treg[i][TIME-1],  m_dot);
	}
	// #endif

/*
		for(i = 0; i < SPACE; i++){
		printf("%g %g %g %g\n",
		Treg[i][n4-1], Treg[i][n2-1], Treg[i][n4], Treg[i][n2-2]);
		}
*/
}



void gnuplot_out(double Treg[SPACE][TIME])
{
	//== gnuplot -> splot =============================
	//-- gnuplot mesh,time mesh�����꤯·����----------
	//-------------------------------------------------
	//-- set dgrid3d [num],[num] - ���ޤ�50,50 --------
	//-- set pm3d -------------------------------------
	//-- splot "file name" u 1:2:3 w d ----------------
	//=================================================
	int i, j;

	for(j = 0; j < TIME-1; j = j + 40){
		for(i = 0; i < SPACE; i = i + 4){
			printf("%d %d %g\n", i, j, Treg[i][j]);
		}
	}
/*
	for(j = 1; j < 501; j += 5){
		for(i = 0; i < 200; i++){
			printf("%d %d %g\n", i, j, Treg[i][j]);
		}
	}
*/
}

int layer_change(int step, int layer_parameter[LAYER+1])
{
	int i;
	int lay;

	for(i = 0; i < LAYER; i++){
//printf("i =%d\n",step);
		if(layer_parameter[i] <= step){
			if(step < layer_parameter[(i+1)]){
				return i;
			}
		}else{
			printf("err layer change!!\n");
			exit(1);
		}
	}
	printf("err layer change2!!\n");
	exit(1);
//	return 0;
}

/********************************************/
/*  imput Temp_C(0T) --> return Capa_C(0T)  */
/*                                          */
/*                                          */
/*          |  T   | C(0T)    | C(5T) |     */
/*          |------|----------|-------|     */
/*          |Temp_A|          |       |     */
/*Temp_C -> |     -|-> Capa_C-|-------|->   */
/*          |Temp_B|          |       |     */
/*          |      |          |       |     */
/*          |      |          |       |     */
/********************************************/

double capacity_old(int column, double data_C[5][15], double imput_temp)
{
	//gas_data[i] = capacity_old(i, He5atm, Tmid)
	int i;
	double capa_A;
	double capa_B;
	double capa_C;
	double temp_A;
	double temp_B;
	double temp_C = imput_temp;

	for(i = 0; i < 15; i++){
		if(temp_C < data_C[0][i]){
			temp_A = data_C[0][i-1];	//1行前の温度データをtemp_Aに代入
			temp_B = data_C[0][i];		//現在(i番目)の温度データをtemp_Bに代入
			capa_A = data_C[column][i-1];	//1行前の(i=1:密度,i=2:比熱,i=3:粘性,i=4:熱伝導係数)データをcapa_Aに代入
			capa_B = data_C[column][i]; //現在(i番目)の温度データをcapa_Bに代入
			break;
		}

		if(i == 14){
			printf("err capacity old!!\n");
			printf("tempC = %f\n", temp_C);
			exit(1);
		}
	}
	capa_C = cal_angle(temp_A, temp_B, temp_C, capa_A, capa_B);
	return capa_C;
}
