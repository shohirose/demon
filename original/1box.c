#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 300//粒子数
#define n ((int)ceil(log2(N)))//完全平衡木の深さ
#define q ((int)pow(2,n)-N)
#define p ((N-q)/2)
double a = 0.02;//radius of disk
double e = 0.95,e_wall = 1.0;//それぞれ粒子同士、壁との反発係数
double g = 1.0;//規格化された重力加速度
double Xmin = -1.0,Xmax = 1.0;//左右の壁の位置
double Ymin = 0.0,Ymax = 1.0;//底面とセルの最高点の位置
double U = 0.149;//床面の振動する速度
double V0 = 1.0;//初期条件での速度分布の標準偏差
int N_cell_x = 32,N_cell_y = 12;//x,y方向のセルの分割数
double T = 20.0;//シミュレーション終了時刻
double epsilon = 0.000001;


struct NODE{//完全平衡木のノードの構造体
	int number;//対応する粒子番号
	double time;//対応する粒子の予測された最短衝突時刻
	struct NODE *left;//上と左右をつなぐ
	struct NODE *right;
	struct NODE *parent;
};

struct EVENT{//ある粒子のイベントの詳細(衝突時刻・相手)を記録
	double time;
	int number_particle;
	int number_col;//今回は使わない
};

struct PARTICLE{//粒子に関する情報をまとめる
	double x,y,u,v;//位置・速度
	double tau;//粒子の固有時間を記録，Delayed State Algorithm(DSA)による高速化のために必要
	double next;//自分の次に同じセルに入った粒子の番号
	struct EVENT event;//次に予定されるイベント
};

struct CELL{
	int first;//セルへの登録時に最も早く登録された粒子
	int last;
};

int intpow(int a,int b);
struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i);
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]);
void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new);
void status_initialize(struct PARTICLE particle[N]);
void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N],int i);
double r_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double v_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x,double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y,double cell_length_y);
void Free_evolution(struct PARTICLE *particle,double t);
void G1(struct PARTICLE *particle,int j);
void G2(struct PARTICLE *particle1,struct PARTICLE *particle2);
double T_DWC(struct PARTICLE particle,double t,int j);
double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t);
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current);
double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max);
double Vmax(struct PARTICLE particle[N]);
void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t);
double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max);


int main(void){
	FILE *fp_position,*fp_height;//ファイルの生成
	char name_position[256];
	sprintf(name_position,"position(N=%d).txt",N);//一定時刻ごとに粒子の位置を保存
	if((fp_position = fopen(name_position,"w"))==NULL){
		printf("file_check open error\n");
	}
	if((fp_height = fopen("height.txt","w"))==NULL){//粒子の高さの平均値を記録
		printf("file open error\n");
	}
	int i_current,j_current;//現在注目している粒子のペア,j_current<0:壁,j_current>=0:粒子
	double v_max = 0.0;//最大速度を保存、セルの更新のために必要
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y =(Ymax-Ymin)/(double)N_cell_y;
	double t=0.0,dt=0.01,trec=0.0,dtrec = (double)T/200.0;//200枚の画像からgifを生成
	double t_cell=0.0,t_cell_old=0.0;//セルの更新時刻
	double height;
	
	//srand((unsigned) time(NULL));
	struct PARTICLE particle[N];
	struct CELL cell[N_cell_x][N_cell_y];
	status_initialize(particle);//位置や速度の初期化
	v_max = Vmax(particle);
	t_cell = (cell_length_y-2.0*a)/(2.0*v_max);//この時間までにマスク外からの衝突はありえない
	cell_register(particle,cell);//粒子をセルに登録する,nextofの初期化
	
	struct NODE *node[n+1][2*p+2*q];//完全平衡木(あるいはトーナメント)を表す構造体
	//nodeの実体化、もうちょっとおしゃれにしたい
	for(int i=0;i<=n;i++){
		for(int j=0;j<2*p+2*q;j++){
			node[i][j] = (struct NODE *)malloc(sizeof(struct NODE));
		}
	}
	
	for(int i=0;i<N;i++){
		particle[i].event = Predictions(particle,cell,t,i);//それぞれの粒子の次のイベントを予測
	}
	CBT_build(node,particle);//Complete Binary Treeを組み立てる
	printf("set up ok\n");
	while(t <= T){
		
		//NEXT EVENTの検索
		i_current = node[0][0]->number;//決勝のノードは最短の時間で衝突する粒子を示す
		j_current = particle[i_current].event.number_particle;//i_currentの衝突相手(これは壁の可能性もある)
		t = NextEvent(particle,cell,node,i_current,j_current);//NEXT EVENTを処理しtとparticle,cell,nodeを更新
		t_cell = t_cell_update(particle[i_current],j_current,t_cell_old,&v_max);//t_cellとv_maxの更新
		
		//i_current,j_currentと衝突する粒子がいた場合はその粒子のeventはinvalidになってしまうので新しくeventを作る
		//そのような粒子は同じマスク内にしか存在しないはずなのでその中で探索
		MaskUpdate(particle,cell,node,i_current,t);//i_currentの周りの粒子でinvalidなものがあればアップデート
		if(j_current >= 0){//jについても同様
			MaskUpdate(particle,cell,node,j_current,t);
		}
		
		
		//EEPGM マスク外の粒子とも衝突する可能性が生じるので登録し直す
		if(t >= t_cell){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			//床に粒子がめり込んでいたらこのエラーが生じる
			for(int i=0;i<N;i++){
				if(particle[i].y < Ymin+a-epsilon){
					printf("i=%d:error\n",i);
					printf("%lf %lf %lf %lf\n",particle[i].x,particle[i].y,particle[i].u,particle[i].v);
					printf("%lf %d %d\n",particle[i].event.time,particle[i].event.number_particle,particle[i].event.number_col);
					G1(&particle[i],-3);
					particle[i].event = Predictions(particle,cell,t,i);
					CBT_update(node,particle[i].event.time,i);
					MaskUpdate(particle,cell,node,i,t);
				}
			}
			
		}
		//粒子の位置の出力
		if((t > trec)&&(t < T)){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			printf("t = %lf, v_max = %lf\n",t,v_max);
			height = 0.0;
			for(int i=0;i<N;i++){
				fprintf(fp_position,"%lf %lf\n",particle[i].x,particle[i].y);
				height += particle[i].y/(double)N;
			}
			fprintf(fp_position,"\n\n");
			fprintf(fp_height,"%lf %lf\n",t,height);
			trec += dtrec;
		}
	}
	
	fclose(fp_position);
	fclose(fp_height);
	return 0;
}

int intpow(int a,int b){
	return (int)pow(a,b);
}

struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i){
	double t_min = 2.0*T,t_temp;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int j_col,j;
	struct PARTICLE particle_j;
	struct EVENT L;
	
	for(j=-3;j<0;j++){
		t_temp = T_DWC(particle[i],particle[i].tau,j);
		if((t_temp > t) & ( t_temp < t_min)){
			t_min = t_temp;
			j_col = j;
		}
	}
	
	int cell_x = getcell_x(particle[i].x,cell_length_x),cell_y = getcell_y(particle[i].y,cell_length_y);
	for(int c1=-1;c1<1;c1++){
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 <= N_cell_x-1)) && (cell_y+c2 >= 0)) && (cell_y+c2 <= N_cell_y-1)){
				j = cell[cell_x+c1][cell_y+c2].first;
				while(j >= 0){
					particle_j = particle[j];
					Free_evolution(&particle_j,particle[i].tau-particle[j].tau);
					t_temp = T_DDC(particle[i],particle[j],particle[i].tau);
					if((t_temp > t) && (t_temp < t_min)){
						t_min = t_temp;
						j_col = j;
					}
					j = particle[j].next;
				}
			}
		}
	}
	L.time = t_min;
	L.number_particle = j_col;
	L.number_col = 0;//今後修正が必要になるかもしれない
	return L;
}
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]){
	int i,n_index;
	//initialization for bottom nodes
	for(i=0;i<2*p+2*q;i++){
		if(i < 2*p){
			node[n][i]->time = particle[i].event.time;
			node[n][i]->number = i;
		}else{
			node[n][i]->time = particle[2*p+(i-2*p)/2].event.time;
			node[n][i]->number = 2*p+(i-2*p)/2;
		}
	}
	//tournament
	for(n_index=n-1;n_index>=0;n_index--){
		for(i=0;i<=intpow(2,n_index)-1;i++){
			node[n_index][i]->left = node[n_index+1][2*i];
			node[n_index+1][2*i]->parent = node[n_index][i];
			node[n_index][i]->right = node[n_index+1][2*i+1];
			node[n_index+1][2*i+1]->parent = node[n_index][i];
			if(node[n_index+1][2*i]->time <= node[n_index+1][2*i+1]->time){
				node[n_index][i]->time = node[n_index+1][2*i]->time;
				node[n_index][i]->number = node[n_index+1][2*i]->number;
			}else{
				node[n_index][i]->time = node[n_index+1][2*i+1]->time;
				node[n_index][i]->number = node[n_index+1][2*i+1]->number;
			}
		}
	}
}

void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new){
	struct NODE *entry_now,hoge_now;
	if(i_new < 2*p){
		entry[n][i_new]->time = time_new;
		entry_now = entry[n][i_new];
	}else{
		entry[n][2*i_new-2*p]->time = time_new;
		entry[n][2*i_new-2*p+1]->time = time_new;
		entry_now = entry[n][2*i_new-2*p];
	}
	while(entry_now->parent != NULL){
		entry_now = entry_now->parent;
		if(entry_now->left->time < entry_now->right->time){
			entry_now->time = entry_now->left->time;
			entry_now->number = entry_now->left->number;
		}else{
			entry_now->time = entry_now->right->time;
			entry_now->number = entry_now->right->number;
		}
	}
}

void status_initialize(struct PARTICLE particle[N]){
	double prob;
	int i;
	for(i=0;i<N;i++){
		prob = Uniform();
		particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
		particle[i].y = (Ymin+a)*prob+(0.5*Ymax-a)*(1-prob);
		while(set(particle,i) == 0){
			prob = Uniform();
			particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+a)*prob+(Ymax-a)*(1-prob);
		}
		particle[i].u = rand_normal(0.0,V0);
		particle[i].v = rand_normal(0.0,V0);
		particle[i].next = -1;
		particle[i].tau = 0.0;
		particle[i].event.number_col = 0;
		particle[i].event.time = 2.0*T;
		particle[i].event.number_particle = -1;
	}
}

void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]){
	int i,cell_x,cell_y,lastPrev;
	//initialize particle.next and cell
	for(i=0;i<N;i++){
		particle[i].next = -1;
	}
	for(cell_x = 0;cell_x < N_cell_x;cell_x++){
		for(cell_y = 0;cell_y < N_cell_y;cell_y++){
			cell[cell_x][cell_y].first = -1;
			cell[cell_x][cell_y].last = -1;
		}
	}

	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	for(i=0;i<N;i++){
		cell_x = getcell_x(particle[i].x,cell_length_x);
		cell_y = getcell_y(particle[i].y,cell_length_y);
		lastPrev = cell[cell_x][cell_y].last;
		cell[cell_x][cell_y].last = i;
		
		if(lastPrev == -1){
			cell[cell_x][cell_y].first = i;
		}else{
			particle[lastPrev].next = i;
		}
	}
}


int set(struct PARTICLE particle[N],int i){//setに成功していれば1,失敗していれば0を返す
	int j,r=1;
	double d;
	
	if(fabs(particle[i].x) < a){
		r = 0;
	}
	for(j=1;j<=i-1;j++){
		d = r_distance(particle[i],particle[j]);
		if(d <= 2.0*a){
			r = 0;
			break;
		}
	}
	return r;
}

double r_distance(struct PARTICLE particle1,struct PARTICLE particle2){
	double d;
	d = sqrt(pow(particle1.x-particle2.x,2.0)+pow(particle1.y-particle2.y,2.0));
	return d;
}

double v_distance(struct PARTICLE particle1,struct PARTICLE particle2){
	double d;
	d = sqrt(pow(particle1.u-particle2.u,2.0)+pow(particle1.v-particle2.v,2.0));
	return d;
}

double Uniform(void){
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}

int getcell_x(double x,double cell_length_x){
	
	if((x < Xmin+a)||(Xmax-a < x)){
		printf("x is out of range\n");
	}
	return (int)((x-Xmin)/cell_length_x);
}

int getcell_y(double y,double cell_length_y){
	if(y < Ymin){
		printf("error:y<0(%lf)\n",y);
		return 0;
	}else if(y>Ymax){
		return N_cell_y-1;//Ymaxよりも高い位置の粒子は一番高いセルに登録
	}else{
		return (int)(y/cell_length_y);
	}
}



void Free_evolution(struct PARTICLE *particle,double t){
	particle->x += (particle->u)*t;
	particle->y += (particle->v)*t-0.5*g*t*t;
	particle->v += -g*t;
	particle->tau += t;//ここはうまくいっているか確認が必要
}

void G1(struct PARTICLE *particle,int j){
	double temp;
	if((j == -1) || (j == -2)){//collision with R or L wall
		particle->u = -e_wall*particle->u;
		if(j == -1){
			particle->x = Xmax-a-epsilon;//このepsilon処理はgetcell_xのときなどに必要になる
		}else{
			particle->x = Xmin+a+epsilon;
		}
	}else if(j == -3){//collision with Bottom wall
		particle->v = (1+e_wall)*U-e_wall*particle->v;
		particle->y = Ymin+a+epsilon;
	}
}


void G2(struct PARTICLE *particle1,struct PARTICLE *particle2){
	double d,Xtemp,Ytemp,Utemp1,Utemp2,Vtemp1,Vtemp2,Cx,Cy;
	d = r_distance(*particle1,*particle2);
	Utemp1 = particle1->u;
	Vtemp1 = particle1->v;
	Utemp2 = particle2->u;
	Vtemp2 = particle2->v;
	Cx = (particle1->x-particle2->x)/d;
	Cy = (particle1->y-particle2->y)/d;
	particle1->u = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
	particle1->v = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
	particle2->u = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
	particle2->v = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
	
}

double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t){
	double r_relative,v_relative,b,hoge;
	double tau = t;
	double x1 = particle1.x ,x2 = particle2.x ,y1 = particle1.y , y2 = particle2.y;
	double u1 = particle1.u ,u2 = particle2.u ,v1 = particle1.v , v2 = particle2.v;
	r_relative = r_distance(particle1,particle2);
	v_relative = v_distance(particle1,particle2);
	b = (x1-x2)*(u1-u2)+(y1-y2)*(v1-v2);
	hoge = b*b-v_relative*v_relative*(r_relative*r_relative-4.0*a*a);
	if(hoge > 0){
		tau += -(b+sqrt(hoge))/(v_relative*v_relative);
	}else{
		tau += T;
	}
	return tau;
}

double T_DWC(struct PARTICLE particle,double t,int j){
	double tau = t;
	if(j==-1){//collision with RIGHT wall(-1)
		if(particle.u>0.0){
			tau += (Xmax-a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-2){//collision with LEFT wall(-2)
		if(particle.u<0.0){
			tau += (Xmin+a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-3){//collision with BOTTOM wall(-3)
		tau += (particle.v+sqrt(particle.v*particle.v+2*g*(particle.y-a)))/g;
	}
	if(tau < t){
		return 2.0*T;
	}else{
		return tau;
	}
}
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current){
	double t = particle[i_current].event.time;
	Free_evolution(&particle[i_current],t-particle[i_current].tau);//i_currentの時間発展
	if(j_current >= 0){//Disk Disk Collision
		Free_evolution(&particle[j_current],t-particle[j_current].tau);//j_currentの時間発展
		G2(&particle[i_current],&particle[j_current]);//粒子同士の衝突処理
	}
	if(j_current < 0){//Disk Wall Collision
		G1(&particle[i_current],j_current);//壁との衝突処理
	}
	particle[i_current].event = Predictions(particle,cell,t,i_current);//i_currentのイベント更新
	CBT_update(node,particle[i_current].event.time,i_current);//i_currentのnodeアップデート
	if(j_current >= 0){//j_currentについても同様
		particle[j_current].event = Predictions(particle,cell,t,j_current);
		CBT_update(node,particle[j_current].event.time,j_current);
	}
	return t;
}

double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max){
	double t_cell,dt_cell;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	if(j_current == -3){
		if(*v_max*(*v_max) < pow(particle.u,2.0)+pow(particle.v,2.0)){
			*v_max = sqrt(pow(particle.u,2.0)+pow(particle.v,2.0));
		}
	}
	dt_cell = (cell_length_y-2.0*a)/(2.0*(*v_max));
	t_cell = t_cell_old+dt_cell;
	return t_cell;
}

double Vmax(struct PARTICLE particle[N]){
	double v_max = 0.0;
	for(int i=0;i<N;i++){
		if(v_max*v_max < particle[i].u*particle[i].u+particle[i].v*particle[i].v){
			v_max = sqrt(particle[i].u*particle[i].u+particle[i].v*particle[i].v);
		}
	}
	return v_max;
}

void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t){
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int cell_x = getcell_x(particle[i_current].x,cell_length_x) , cell_y = getcell_y(particle[i_current].y,cell_length_y),j;
	for(int c1=-1;c1<=1;c1++){
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)){
				j = cell[cell_x+c1][cell_y+c2].first;
				while(j >= 0){
					if(particle[j].event.number_particle == i_current){
						particle[j].event = Predictions(particle,cell,t,j);
						CBT_update(node,particle[j].event.time,j);
					}
					j = particle[j].next;
				}
			}
		}
	}
}

double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max){
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	double dt_cell,t_cell;

	for(int i=0;i<N;i++){//現在時刻まで時間発展
		Free_evolution(&particle[i],t-particle[i].tau);
	}

	cell_register(particle,cell);//全粒子をセルに登録し直す
	for(int i=0;i<N;i++){//全粒子についてeventを計算し直す
		particle[i].event = Predictions(particle,cell,t,i);
	}
	CBT_build(node,particle);//CBTも最初から構成
	*v_max = Vmax(particle);
	dt_cell = (cell_length_y-2.0*a)/(2.0*(*v_max));
	t_cell = t+dt_cell;
	return t_cell;
}