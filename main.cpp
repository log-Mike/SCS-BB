//  
//  main.cpp
//  CSSC
//
//  Created by kai on 25/5/2020.
//  Copyright © 2020 kai. All rights reserved.
//

#include "Utility.h"
#include "ListLinearHeap.h"
#include "Timer.h"
#include <fstream>

ofstream fout;

int MaxTime = 1800;

ui domBr;
ui binBr;

double StartTime = 0.0;

bool EXE_heu2;
bool EXE_heu3;
bool EXE_heu4;

bool EXE_ub1;
bool EXE_ub2;
bool EXE_ub3;
bool EXE_ub3_optimization;
bool EXE_refine_G0;
bool EXE_core_maintenance;
bool EXE_new2VI;
bool EXE_del_from_VR;
bool EXE_dom_ustar;
double total_val_ub1;
double total_val_ub3;
double total_UB;
ui domS_Threshold;
ui srch_ord;
double total_Heu_time;
bool over_time_flag;

ui n; //vertices
ui m; //edges
ui l; //size LB
ui h; //size UB
int QID; // vertex we need in solution
ui dMAX; // maximum degree
ui kl; //min min deg
ui ku; //max min deg
vector<ui> H; // current optimal solution
ui ubD = INF; // upper bound on optimal diameter

ui * pstart;
ui * edges;
ui * peel_sequence;
ui * degree;
ui * core;
ui * q_dist;

vector<ui> G0;
ui * G0_edges;
ui * G0_x;
ui * G0_deg;

vector<ui> VI;
vector<ui> VIVR;
bool * inVI;
bool * inVR;
ui * degVI;
ui * degVIVR;

vector<ui> NEI;
ui * inNEI;
double * NEI_score;
vector<vector<ui>> combs;

/*
defined never used 
bool cmp_of_domS(const ui x, const ui y)
{
    return degVIVR[x]>degVIVR[y];
}*/

double time_new2VI;
double time_del_from_VR;
double time_find_NEI;
double time_find_usatr;
double time_comp_ub;

void load_graph(const char * graph_file)
{
    map<ui, set<ui>> G;
    string buffer;
    ifstream inputFile(graph_file, ios::in);
    
    // load graph as AL into G
    if (!inputFile.is_open())
    {
        cout << "Graph file Open Failed "<<endl;
        exit(1);
    }
    else
    {
        inputFile >> n >> m;
        
        ui tu,tv;
        while (inputFile >> tu >> tv)
        {
            if(tu==tv) continue;
            G[tu].insert(tv);
            G[tv].insert(tu);
        }
        inputFile.close();
    }
    

    // set peel sequence to i for
    // computing coreness later
    // set node degrees
    // find max degree
    peel_sequence = new ui[n];
    degree = new ui[n];
    core = new ui[n];
    dMAX = 0;

    for(ui i = 0; i<n; i++){
        peel_sequence[i] = i;
        // empty AL -> end -> deg 0
        if(G.find(i)!=G.end()){
            degree[i] = G.find(i)->second.size();
            if(degree[i]>dMAX) dMAX = degree[i];
        }
        else degree[i] = 0;
    }



    // storing graph as 2 
    // 1d arrays as well
    pstart = new ui[n+1];
    edges = new ui[2*m];

    pstart[0] = 0;
    for(ui i =0;i<n;i++){
        if(G.find(i)!=G.end()){
            ui j = 0;
            for(ui nei : G[i]){
                edges[pstart[i]+j] = nei;j++;
            }
            pstart[i+1] = pstart[i] + G[i].size();
        }
        else{
            pstart[i+1] = pstart[i];
        }
    }
    cout<<"n="<<n<<",m="<<m<<",dMAX="<<dMAX<<endl;
}

FILE *open_file(const char *file_name, const char *mode) {
    FILE *f = fopen(file_name, mode);
    if(f == nullptr) {
        printf("Can not open file: %s\n", file_name);
        exit(1);
    }

    return f;
}

/*
defined never used
void load_graph_binv(const char* dir)
{
    printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
    FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

    ui tt;
    fread(&tt, sizeof(ui), 1, f);
    if(tt != sizeof(ui)) {
        printf("sizeof unsigned int is different: b_degree.bin(%u), machine(%lu)\n", tt, sizeof(ui));
        return ;
    }

    fread(&n, sizeof(ui), 1, f);
    fread(&m, sizeof(ui), 1, f);
    m=m/2;
    cout<<"binv, n="<<n<<", m="<<m<<endl;

    degree = new ui[n];//degree
    
    fread(degree, sizeof(ui), n, f);//size_t fread ( void *buffer, size_t size, size_t count, FILE *stream) ;

    fclose(f);

    f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

    if(pstart != nullptr) delete[] pstart;
    pstart = new ui[n+1];
    if(edges != nullptr) delete[] edges;
    edges = new ui[2*m];
    
    peel_sequence = new ui[n];
    core = new ui[n];
    dMAX = 0;
    for(ui i = 0; i<n; i++){
        peel_sequence[i] = i;
        if(degree[i]>dMAX) dMAX = degree[i];
    }
    pstart[0] = 0;
    for(ui i = 0;i < n;i ++) {
        if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
        pstart[i+1] = pstart[i] + degree[i];
    }

    fclose(f);
}*/

/*void GreedyDist()
defined never used
{
    bool * bit = new bool[n];
    memset(bit, 0, sizeof(bool)*n);
    bit[QID] = 1;
    for(ui i=pstart[QID]; i < pstart[QID+1]; i++){
        ui nei = edges[i];
        bit[nei] = 1;
    }
    ui mindeg = degree[QID];
    for(ui i = pstart[QID]; i < pstart[QID+1]; i++){
        ui nei = edges[i];
        ui v_d = 0;
        for(ui j = pstart[nei]; j < pstart[nei+1]; j++){
            ui v = edges[j];
            if(bit[v])
                ++ v_d;
        }
        if(v_d < mindeg)
            mindeg = v_d;
    }
    delete [] bit;
}*/

void core_decomposition_linear_list()
{
    ui max_core = 0;
    ListLinearHeap *linear_heap = new ListLinearHeap(n, n-1);
    linear_heap->init(n, n-1, peel_sequence, degree);
    memset(core, 0, sizeof(ui)*n);
    for(ui i = 0;i < n;i ++) {
        ui u, key;
        linear_heap->pop_min(u, key);
        if(key > max_core) max_core = key;
        peel_sequence[i] = u;
        core[u] = max_core;
        for(ui j = pstart[u];j < pstart[u+1];j ++) if(core[edges[j]] == 0) {
            linear_heap->decrement(edges[j]);
        }
    }
    delete linear_heap;
}

void cal_query_dist()
{
    q_dist = new ui[n];
    for(ui i =0;i<n;i++)
        q_dist[i] = INF;
    queue<ui> Q;
    q_dist[QID] = 0;
    Q.push(QID);
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        for(ui i = pstart[v]; i < pstart[v+1]; i++){
            ui w = edges[i];
            if(q_dist[w] == INF){
                q_dist[w] = q_dist[v] + 1;
                Q.push(w);
            }
        }
    }
}

/*void heu1(vector<ui> & H1, ui & kl1)
defined never used
{
    H1.clear();
    kl1 = 0;
    
    vector<ui> tH;
    tH.clear();
    ui tkl = 0;
    ui tH_size=0;
    
    ui * sta = new ui[n];
    memset(sta, 0, sizeof(ui)*n);
    
    ui * deg = new ui[n];
    memset(deg, 0, sizeof(ui)*n);
    
    priority_queue<pair<ui, ui>> Q;
    Q.push(make_pair(degree[QID], QID));
    sta[QID] = 1;
    while (!Q.empty()) {
        ui v = Q.top().second;
        Q.pop();
        sta[v] = 2;
        tH.push_back(v);
        for(ui nei = pstart[v]; nei<pstart[v+1]; nei++){
            if(!sta[edges[nei]]){
                Q.push(make_pair(degree[edges[nei]], edges[nei]));
                sta[edges[nei]] = 1;
            }
            if(sta[edges[nei]]==2){
                ++ deg[edges[nei]];
                ++ deg[v];
            }
        }
        if(tH.size()>=l ){
            ui cur_mind = INF;
            for(ui i=0; i<tH.size(); i++){
                if(deg[tH[i]]<cur_mind)
                    cur_mind = deg[tH[i]];
            }
            if(cur_mind>=tkl){
                tkl = cur_mind;
                tH_size = tH.size();
            }
        }
        if(tH.size()==h) break;
    }
    
    for(ui i=0; i<tH_size; i++){
        H1.push_back(tH[i]);
    }
    kl1 = tkl;
    
    delete [] sta;
    delete [] deg;
}*/

void heu2(vector<ui> & H2, ui & kl2)
{
    H2.clear();
    kl2 = 0;
    
    vector<ui> tH;
    ui tkl = 0;
    ui tH_size = 0;
    
    ui * sta = new ui[n];
    memset(sta, 0, sizeof(ui)*n);
    
    ui * deg = new ui[n];
    memset(deg, 0, sizeof(ui)*n);
    
    priority_queue<pair<ui, ui>> Q;
    Q.push(make_pair(0, QID));
    sta[QID] = 1;
    
    while(!Q.empty()){
        ui v = Q.top().second;
        Q.pop();
        if(sta[v] == 2) continue;
        tH.push_back(v);
        sta[v] = 2;
        for(ui j = pstart[v]; j < pstart[v+1]; j++){
            ui nei = edges[j];
            if(sta[nei] == 0){
                ui d = 0;
                for(ui w = pstart[nei]; w < pstart[nei+1]; w++){
                    if(sta[edges[w]] == 2) ++ d;
                }
                Q.push(make_pair(d, nei));
                sta[nei] = 1;
            }
            else{
                if(sta[nei] == 1){
                    ui new_d = 0;
                    for(ui w = pstart[nei]; w < pstart[nei+1]; w++){
                        if(sta[edges[w]] == 2) ++ new_d;
                    }
                    Q.push(make_pair(new_d, nei));
                }
                else{
                    ++ deg[nei];
                    ++ deg[v];
                }
            }
        }
        if(tH.size()>=l){
            ui mindeg = INF;
            for(ui i = 0; i < tH.size(); i++){
                if(deg[tH[i]] < mindeg)
                    mindeg = deg[tH[i]];
            }
            if(mindeg >= tkl){
                tkl = mindeg;
                tH_size = (ui)tH.size();
            }
        }
        if(tH.size()==h)
            break;
    }
    
    for(ui i=0; i<tH_size; i++){
        H2.push_back(tH[i]);
    }
    kl2 = tkl;
    
    delete [] sta;
    delete [] deg;
}

void heu3(vector<ui> & H3, ui & kl3)
{
    H3.clear();
    kl3 = 0;
    
    vector<ui> tH;
    ui tkl = 0;
    ui tH_size = 0;
    
    ui * sta = new ui[n];
    memset(sta, 0, sizeof(ui)*n);
    
    ui * deg = new ui[n];
    memset(deg, 0, sizeof(ui)*n);
    
    priority_queue<pair<double, ui>> Q;
    Q.push(make_pair(0, QID));
    sta[QID] = 1;
    
    while (!Q.empty()) {
        ui v = Q.top().second;
        Q.pop();
        if(sta[v] == 2) continue;
        tH.push_back(v);
        sta[v] = 2;
        for(ui nei = pstart[v]; nei<pstart[v+1]; nei++){
            if(sta[edges[nei]] == 2){
                ++ deg[edges[nei]];
                ++ deg[v];
            }
        }
        for(ui nei = pstart[v]; nei<pstart[v+1]; nei++){
            if(!sta[edges[nei]]){
                double score = 0;
                for(ui w = pstart[edges[nei]]; w < pstart[edges[nei]+1]; w++){
                    if(sta[edges[w]] == 2 && deg[edges[w]] != 0){
                        score += (double) 1/deg[edges[w]];
                    }
                }
                score += (double) degree[edges[nei]]/dMAX;
                Q.push(make_pair(score, edges[nei]));
                sta[edges[nei]] = 1;
            }
            else{
                if(sta[edges[nei]] == 1){
                    double new_score = 0;
                    for(ui w = pstart[edges[nei]]; w < pstart[edges[nei]+1]; w++){
                        if(sta[edges[w]] == 2 && deg[edges[w]] != 0){
                            new_score += (double) 1/deg[edges[w]];
                        }
                    }
                    new_score += (double) degree[edges[nei]]/dMAX;
                    Q.push(make_pair(new_score, edges[nei]));
                }
            }
        }
        if(tH.size()>=l){
            ui mindeg = INF;
            for(ui i = 0; i < tH.size(); i++){
                if(deg[tH[i]] < mindeg)
                    mindeg = deg[tH[i]];
            }
            if(mindeg >= tkl){
                tkl = mindeg;
                tH_size = (ui)tH.size();
            }
        }
        if(tH.size()==h)
            break;
    }
    
    for(ui i=0; i<tH_size; i++){
        H3.push_back(tH[i]);
    }
    kl3 = tkl;
    
    delete [] sta;
    delete [] deg;
}

void heu4(vector<ui> & H4, ui & kl4)
{
    H4.clear();
    kl4 = 0;
    
    if( degree[QID] < l-1){
        return;
    }
    
    vector<ui> S;
    bool * inS = new bool[n];
    memset(inS, 0, sizeof(bool)*n);
    
    ui * deg = new ui[n];
    memset(deg, 0, sizeof(ui)*n);
    
    S.push_back(QID);
    inS[QID] = 1;
    ++ deg[QID];
    
    for(ui j = pstart[QID]; j < pstart[QID+1]; j++){
        S.push_back(edges[j]);
        inS[edges[j]] = 1;
        for(ui w = pstart[edges[j]]; w < pstart[edges[j]+1]; w++){
            if(inS[edges[w]]){
                ++ deg[edges[w]];
                ++ deg[edges[j]];
            }
        }
    }

    priority_queue<pair<ui,ui>, vector<pair<ui,ui>>, greater<>> Q;
    for(auto e : S)
        Q.push(make_pair(deg[e], e));
    vector<ui> rmv;
    ui rmv_idx = 0;
    ui mindeg = 0;
    while (!Q.empty()) {
        ui d = Q.top().first;
        ui v = Q.top().second;
        Q.pop();
        if(inS[v] == 0) continue;
        ui remain = (ui)S.size() - (ui)rmv.size();
        if(remain >= l && remain <= h){
            if(d > mindeg){
                mindeg = d;
                rmv_idx = (ui)rmv.size();
            }
        }
        if(remain == l) break;

        inS[v] = 0;
        rmv.push_back(v);
        for(ui j = pstart[v]; j < pstart[v+1]; j++){
            if(inS[edges[j]]){
                -- deg[edges[j]];
                Q.push(make_pair(deg[edges[j]], edges[j]));
            }
        }
    }
    
    memset(inS, 0, sizeof(bool)*n);
    for(auto e : S)
        inS[e] = 1;
    for(ui i = 0; i < rmv_idx; i++)
        inS[rmv[i]] = 0;
    
    vector<ui> tH;
    ui tkl = mindeg;
    for(auto e : S)
        if(inS[e]) tH.push_back(e);
    
    for(ui i=0; i<tH.size(); i++){
        H4.push_back(tH[i]);
    }
    kl4 = tkl;
    
    delete [] inS;
    delete [] deg;
    
}

/*
defined never used
void heu5(vector<ui>& H5, ui& kl5)
{
    H5.clear();
    kl5 = 0;
    
    vector<ui> S;
    bool * inS = new bool[n];
    memset(inS, 0, sizeof(bool)*n);
    
    ui * deg = new ui[n];
    memset(deg, 0, sizeof(ui)*n);
    
    inS[QID] = 1;
    S.push_back(QID);
    for(ui j = pstart[QID]; j <pstart[QID+1]; j++){
        if(!inS[edges[j]]){
            inS[edges[j]] = 1;
            S.push_back(edges[j]);
        }
        for(ui w = pstart[edges[j]]; w < pstart[edges[j]+1]; w++){
            if(!inS[edges[w]]){
                inS[edges[w]] = 1;
                S.push_back(edges[w]);
            }
        }
    }
    
    if(S.size() < l){
        return;
    }
    
    for(auto e : S){
        for(ui i = pstart[e]; i < pstart[e+1]; i++){
            if(inS[edges[i]]) ++deg[edges[i]];
        }
    }
    
    priority_queue<pair<ui,ui>, vector<pair<ui,ui>>, greater<>> Q;
    for(auto e : S)
        Q.push(make_pair(deg[e], e));
    vector<ui> rmv;
    ui rmv_idx = 0;
    ui mindeg = 0;
    while (!Q.empty()) {
        ui v = Q.top().second;
        ui d = Q.top().first;
        Q.pop();
        if(!inS[v] || v==QID) continue;

        ui remain = (ui)S.size() - (ui)rmv.size();
        if(remain >= l && remain <= h){
            ui rd = min(deg[QID], d);
            if(rd >= mindeg){
                mindeg = rd;
                rmv_idx = (ui)rmv.size();
            }
        }
        if(remain == l) break;
        inS[v] = 0;
        rmv.push_back(v);
        for(ui i = pstart[v]; i < pstart[v+1]; i++){
            if(inS[edges[i]]){
                -- deg[edges[i]];
                Q.push(make_pair(deg[edges[i]], edges[i]));
            }
        }
    }
    memset(inS, 0, sizeof(bool)*n);
    for(auto e : S)
        inS[e] = 1;
    for(ui i = 0; i < rmv_idx; i++)
        inS[rmv[i]] = 0;
    
    vector<ui> tH;
    ui tkl = mindeg;
    for(auto e : S)
        if(inS[e]) tH.push_back(e);
    
    for(ui i=0; i<tH.size(); i++){
        H5.push_back(tH[i]);
    }
    kl5 = tkl;
    
    delete [] inS;
    delete [] deg;
}*/

void CSSC_heu()
{
    H.clear();
    kl = 0;
    
    vector<ui> H2; ui kl2 = 0;
    if(EXE_heu2) heu2(H2,kl2);
    vector<ui> H3; ui kl3 = 0;
    if(EXE_heu3) heu3(H3,kl3);
    vector<ui> H4; ui kl4 = 0;
    if(EXE_heu4) heu4(H4,kl4);
    
    if(kl4 >= kl3 && kl4 >= kl2){
        H = H4;
        kl = kl4;
    }
    else{
        if(kl3 >= kl2){
            H = H3;
            kl = kl3;
        }
        else{
            H = H2;
            kl = kl2;
        }
    }
}

void reduction_g()
{
    G0.clear();
    G0_edges = new ui[2*m];
    
    // used interchangably
    // for degree
    G0_x = new ui[n];
    G0_deg  = new ui[n];
    
    memset(G0_x, 0, sizeof(ui)*n);
    memset(G0_deg, 0, sizeof(ui)*n);
    
    bool * inQ = new bool[n];
    memset(inQ, 0, sizeof(bool)*n);
    
    queue<ui> Q;
    Q.push(QID);
    inQ[QID] = 1;
    
    // BFS from query vertex
    // only add nodes not removable 
    // by cn(v) <= kl

   
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        G0.push_back(v);
        for(ui i = pstart[v]; i < pstart[v+1]; i++){
            // doesn't get to lower bound or is out of distance:
            // (if distance from query is farther than upper bound
            // on diameter, can safely remove it)
            if(core[edges[i]] > kl && q_dist[edges[i]] <= ubD){
                G0_edges[pstart[v] + G0_x[v]] = edges[i];
                ++ G0_x[v];
                ++ G0_deg[v];
                if(!inQ[edges[i]]){
                    Q.push(edges[i]);
                    inQ[edges[i]] = 1;
                }
            }
        }
    }
    // resulting graph:
    // p_start (starting index of edges for a vertex) are the same
    // G0_edges is the same length but will contain much less vertices
    // this is fine since we know degree

    // instead of:
    // i = pstart[v]; i < pstart[v + 1]
    // it is now:
    // i = pstart[v]; i < pstart[v] + G0_x[v]
    delete [] inQ;
}

/*
defined never used
void refine_G0()
{
//    cout<<"in refine_G0"<<endl;
    queue<ui> rmv;
    bool * inrmv = new bool[n];
    memset(inrmv, 0, sizeof(bool)*n);
    
    for(auto e : VIVR){
        if(degVIVR[e] <= kl){
            rmv.push(e);
            inrmv[e] = 1;
        }
    }
    
    while (!rmv.empty()) {
        ui v = rmv.front();
        rmv.pop();
        inVR[v] = 0;
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVR[w]){
                -- degVIVR[w];
                -- degVIVR[v];
//                cout<<" -- degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//                cout<<" -- degVIVR["<<v<<"]="<<degVIVR[v]<<endl;
                if(degVIVR[w] <= kl && !inrmv[w]){
                    rmv.push(w);
                    inrmv[w] = 1;
                }
            }
        }
    }
    delete [] inrmv;
}*/

// find u*
int find_ustar()
{
    int uid = -1;
    double best_score = 0;
    // NEI stores 1 hop neighbors in candidate set
    // pull from there for u*
    for(auto e : NEI){
        if(inVR[e]){
            double its_score = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if(inVI[G0_edges[i]] && degVI[G0_edges[i]] != 0)
                  its_score += (double) 1/degVI[G0_edges[i]];
            }
            its_score += (double)degVIVR[e]/dMAX;
            if(its_score > best_score){
                best_score = its_score;
                uid = e;
            }
        }
    }
    return uid;
}

// same as regular find u* but start search with neighbors
// of lowest degree vertices in partial solution
int find_ustar_mindeg()
{
    double best_score;
    int uid = -1;
    set<ui> dict_deg;
    for(auto e : VI){
        dict_deg.insert(degVI[e]);
    }

    for(auto deg : dict_deg){
        // vt stores top 1, 2, etc lowest deg
        // vertices in partial solution
        vector<ui> vt;
        for(auto e : VI){
            if(degVI[e]==deg){
                vt.push_back(e);
            }
        }

        best_score = 0;
        for(auto e : vt){
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i];
                if(inVR[w]){
                    double its_score = 0;
                    for(ui j = pstart[w]; j < pstart[w]+G0_x[w]; j++){
                        if(inVI[G0_edges[j]] && degVI[G0_edges[j]] != 0){
                            its_score += (double) 1/degVI[G0_edges[j]];
                        }
                    }
                    its_score += (double)degVIVR[w]/dMAX;
                    if(its_score > best_score){
                        best_score = its_score;
                        uid = w;
                    }
                }
            }
        }
        if(uid != -1)
            break;
    }
    return uid;
}

// depending on threshhold, find u* with one of two above functions
int find_ustar_2phase()
{
    if( VI.size() > (h*2)/5 )
        return find_ustar_mindeg();
    else
        return find_ustar();
}

// same code as original find u*, not finished ?!?!
int find_ustar_link()
{
    int uid = -1;
    double best_score = 0;
    for(auto e : NEI){
        if(inVR[e]){
            double its_score = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if(inVI[G0_edges[i]])
                  its_score ++;
            }
            if(its_score >= best_score){
                best_score = its_score;
                uid = e;
            }
        }
    }
    return uid;
}

// random vertex as u*, essentially baseline enum
int find_ustar_random()
{
    int uid = -1;
    for(auto e : NEI){
        if(inVR[e]){
            return  e;
        }
    }
    return uid;
}

// finding a min
ui get_ub1()
{
    ui min_deg = INF;
    ui r = h - (ui)VI.size();
    if(r<0){
        cout<<"??? r < 0 ???"<<endl;
        exit(1);
    }
    // find min
    for(auto e : VI){
        ui its_deg_ub = 0;
        ui cands = 0;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(inVR[G0_edges[i]]) ++ cands;
        }
        its_deg_ub = degVI[e] + min(r, cands);
        if(its_deg_ub < min_deg) min_deg = its_deg_ub;
    }
    return min_deg;
}

// top k largest/smallest
ui get_ub2()
{
    vector<ui> nei;
    bool * innei = new bool[n];
    memset(innei, 0, sizeof(bool)*n);
    for(auto e : VI){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(inVR[G0_edges[i]] && !innei[G0_edges[i]]){
                nei.push_back(G0_edges[i]);
                innei[G0_edges[i]] = 1;
            }
        }
    }
    delete [] innei;
    
    ui r = h - (ui)VI.size();
    vector<ui> cov_power;
    for(auto e : nei){
        ui its_power = 0;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(inVI[G0_edges[i]])
                ++ its_power;
        }
        cov_power.push_back(its_power);
    }
    ui r1 = (ui)cov_power.size();
    ui r2 = min(r, r1);
    sort(cov_power.begin(), cov_power.end(), greater<>());//decreasing order
    vector<ui> interm_deg;
    for(auto e : VI){
        interm_deg.push_back(degVI[e]);
    }
    sort(interm_deg.begin(), interm_deg.end(), less<>());//increasing order
    
    for(ui i = 0; i < r2; i++){
        ui budget = cov_power[i];
        for(ui j = 0; j < budget; j++){
            ++ interm_deg[j];
        }
        sort(interm_deg.begin(), interm_deg.end(), less<>());//increasing order
    }
    
    return interm_deg[0];
}

ui get_ub3()
{
    Timer t;
    vector<ui> nei;
    bool * innei = new bool[n];
    memset(innei, 0, sizeof(bool)*n);
    for(auto e : VI){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            if(inVR[G0_edges[i]] && !innei[G0_edges[i]]){
                nei.push_back(G0_edges[i]);
                innei[G0_edges[i]] = 1;
            }
        }
    }
    delete [] innei;
    
    set<ui> ditc_deg;
    for(auto e : VI){
        ditc_deg.insert(degVI[e]);
    }
    vector<ui> deg_lev;
    for(auto e : ditc_deg)
        deg_lev.push_back(e);
    
    ui ditc_deg_num = (ui)deg_lev.size();
    
    ui r = h - (ui)VI.size();
    
    ui min_deg = INF;
    ui rcd_i = 0;
    
    for(ui i = 0; i < ditc_deg_num; i++){
        ui t_deg = deg_lev[i];
        vector<ui> interm_deg;
        for(auto e : VI){
            if(degVI[e] <= t_deg){
                interm_deg.push_back(degVI[e]);
            }
        }
        sort(interm_deg.begin(), interm_deg.end(), less<>()); //increasing order
        
        vector<ui> cov_power;
        for(auto e : nei){
            ui its_power = 0;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                if(inVI[G0_edges[i]] && degVI[G0_edges[i]] <= t_deg){
                    ++ its_power;
                }
            }
            if(its_power > 0){
                cov_power.push_back(its_power);
            }
        }
        sort(cov_power.begin(), cov_power.end(), greater<>()); //decreasing order
        
        ui r1 = (ui)cov_power.size();
        ui r2 = min(r, r1);
        
        for(ui i = 0; i < r2; i++){
            ui budget = cov_power[i];
            for(ui j = 0; j < budget; j++){
                ++ interm_deg[j];
            }
            sort(interm_deg.begin(), interm_deg.end(), less<>());//increasing order
        }
        if(min_deg > interm_deg[0]){
            min_deg = interm_deg[0];
            rcd_i = i;
        }
        if(min_deg <= kl)
            return min_deg;
        
        if(EXE_ub3_optimization){
            if(i < ditc_deg_num-1 && min_deg <= deg_lev[i+1]){
                return min_deg;
            }
        }
    }
    return min_deg;
}

/*ui get_ub33()               Same as ub3 w/clean up on coding
{
    set<ui> ditc_deg;
    for(auto e : VI){
        ditc_deg.insert(degVI[e]);
    }
    vector<ui> deg_lev;
    for(auto e : ditc_deg)
        deg_lev.push_back(e);
    
    ui ditc_deg_num = (ui)deg_lev.size();
    
    ui r = h - (ui)VI.size();
    
    ui min_deg = INF;
    ui rcd_i = 0;
    
    for(ui i = 0; i < ditc_deg_num; i++){
        ui t_deg = deg_lev[i];
//        ui min_deg = INF;
        vector<ui> interm_deg;
        for(auto e : VI){
            if(degVI[e] <= t_deg){
                interm_deg.push_back(degVI[e]);
            }
        }
        sort(interm_deg.begin(), interm_deg.end(), less<>()); //increasing order
        
        vector<ui> cov_power;
        for(auto e : NEI){
            if(inNEI[e] != 0){
                ui its_power = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]] && degVI[G0_edges[i]] <= t_deg){
                        ++ its_power;
                    }
                }
                if(its_power > 0){
                    cov_power.push_back(its_power);
                }
            }
        }
        sort(cov_power.begin(), cov_power.end(), greater<>()); //decreasing order
        
        ui r1 = (ui)cov_power.size();
        ui r2 = min(r, r1);
        
        for(ui i = 0; i < r2; i++){
            ui budget = cov_power[i];
            for(ui j = 0; j < budget; j++){
                ++ interm_deg[j];
            }
            sort(interm_deg.begin(), interm_deg.end(), less<>());//increasing order
        }
        if(min_deg > interm_deg[0]){
            min_deg = interm_deg[0];
            rcd_i = i;
        }
        if(EXE_ub3_optimization){
            if(i < ditc_deg_num-1 && min_deg <= deg_lev[i+1]){
                return min_deg;
            }
        }
    }
    return min_deg;
}*/


// compute upper bound
// takes min(Ub1, Ub2, Ub3)
// Not the same as how they
// described in their paper
ui compute_ub()
{
    ++ total_UB;
    ui ub1 = INF;
    if(EXE_ub1) ub1 = get_ub1();
    ui ub2 = INF;
    if(EXE_ub2) ub2 = get_ub2();
    ui ub3 = INF;
    if(EXE_ub3) ub3 = get_ub3();
    // total_val_ub1 += ub1;        Not used in rest of code
    // total_val_ub3 += ub3;		Not used in rest of code
    return min(ub1, min(ub2,ub3));
}



/*
defined not used
void upd_based_kl(bool & del_v_in_VI, vector<ui> & rVIVR, vector<bool> & rVIVRbit, ui ustar)
{
    inVR[ustar] = 1;
    bool * inQ = new bool[n];
    memset(inQ, 0, sizeof(bool)*n);
    queue<ui> Q;
    Q.push(ustar);
    inQ[ustar] = 1;
    
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        rVIVR.push_back(v);
        if(inVI[v]) {del_v_in_VI = true; inVI[v]=0; rVIVRbit.push_back(1);}
        if(inVR[v]) {inVR[v]=0; rVIVRbit.push_back(0);}
        
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVI[w]){
                --degVIVR[w];
                if(degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);inQ[w] = 1;
                }
            }
            if(inVR[w]){
                --degVIVR[w];
                if(degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);inQ[w] = 1;
                }
            }
        }
    }
    inVR[ustar] = 0;
    delete [] inQ;
}*/

// 
void core_maintenance(bool & del_v_in_VI, vector<ui> & rVR, ui ustar)
{
    
    queue<ui> Q;
    bool * inQ = new bool[n];
    memset(inQ, 0, sizeof(bool)*n);
    
    
    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
        ui w = G0_edges[i];
        if((inVI[w] || inVR[w]) && degVIVR[w] <= kl){
            Q.push(w);
            inQ[w] = 1;
        }
    }
    
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(inVI[v]){
            del_v_in_VI = true;
            break;
        }
        else if(inVR[v]){
            inVR[v] = 0;
            rVR.push_back(v);
        }
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVI[w] || inVR[w]){
                -- degVIVR[w];
                -- degVIVR[v];
//                cout<<" due to the removal of vertex "<<v<<", --degVIVR["<<w<<"] = "<<degVIVR[w]<<endl;
//                cout<<" due to the removal of vertex "<<v<<", --degVIVR["<<v<<"] = "<<degVIVR[v]<<endl;
                if(degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);
                    inQ[w] = 1;
                }
            }
        }
    }
    delete [] inQ;
}

void core_maintenance(bool & del_v_in_VI, vector<ui> & rVR, vector<ui> & del_vec)
{
    queue<ui> Q;
    bool * inQ = new bool[n];
    memset(inQ, 0, sizeof(bool)*n);
    
    for(auto e : del_vec){
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            ui w = G0_edges[i];
            if( (inVI[w] || inVR[w]) && degVIVR[w] <= kl && !inQ[w]){//?????????
                Q.push(w);
                inQ[w] = 1;
            }
        }
    }

    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(inVI[v]){
            del_v_in_VI = true;
            break;
        }
        else if(inVR[v]){
            inVR[v] = 0;
            rVR.push_back(v);
        }
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVI[w] || inVR[w]){
                -- degVIVR[w];
                -- degVIVR[v];
                if(degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);
                    inQ[w] = 1;
                }
            }
        }
    }//while
    delete [] inQ;
}

void core_maintenance(bool & del_v_in_VI, vector<ui> & rVR, vector<pair<double, ui>> & del_vec)
{
    queue<ui> Q;
    bool * inQ = new bool[n];
    memset(inQ, 0, sizeof(bool)*n);
    
    for(auto ee : del_vec){
        auto e = ee.second;
        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
            ui w = G0_edges[i];
            if( (inVI[w] || inVR[w]) && degVIVR[w] <= kl && !inQ[w]){//?????????
                Q.push(w);
                inQ[w] = 1;
            }
        }
    }

    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(inVI[v]){
            del_v_in_VI = true;
            break;
        }
        else if(inVR[v]){
            inVR[v] = 0;
            rVR.push_back(v);
        }
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVI[w] || inVR[w]){
                -- degVIVR[w];
                -- degVIVR[v];
                if(degVIVR[w] <= kl && !inQ[w]){
                    Q.push(w);
                    inQ[w] = 1;
                }
            }
        }
    }//while
    delete [] inQ;
}

/*
defined never used
void upd_based_ubD(bool& del_v_in_VI_ubD, ui ustar)
{
    int * dist = new int[n];
    memset(dist, -1, sizeof(int)*n);
    queue<ui> Q;
//    bool * inQ = new bool[n];
//    memset(inQ, 0, sizeof(bool)*n);
    
    dist[ustar] = 0;
    Q.push(ustar);
//    inQ[ustar] = 1;
    ui max_d = 0;
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(dist[v] > max_d)
            max_d = dist[v];
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if(inVI[w] && dist[w] == -1){
                dist[w] = dist[v] + 1;
                Q.push(w);
            }
        }
    }
    delete [] dist;
    if(max_d >= ubD)
        del_v_in_VI_ubD = true;
}*/

/* void upd_based_ubD1(bool & del_v_in_VI_ubD, ui ustar)
defined never used again
{
    int * dist = new int[n];
    memset(dist, -1, sizeof(int)*n);
    
    bool * be_touched = new bool[n];
    memset(be_touched, 0, sizeof(bool)*n);
    
    queue<ui> Q;
    dist[ustar] = 0;
    Q.push(ustar);
    
    while (!Q.empty()) {
        ui v = Q.front();
        Q.pop();
        if(dist[v] > ubD){
            break;
        }
        be_touched[v] = 1;
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            ui w = G0_edges[i];
            if( (inVR[w] || inVI[w]) && dist[w] == -1 ){
                dist[w] = dist[v] + 1;
                Q.push(w);
            }
        }
    }
    
    for(auto e : VI){
        if(!be_touched[e]){
            del_v_in_VI_ubD = true;
            break;
        }
    }
    delete [] dist;
    delete [] be_touched;
}*/

/*
defined never used
void upd_based_ubD(bool& del_v_in_VI_ubD, vector<ui>& add_vec)
{
    int * dist = new int[n];
    queue<ui> Q;
    ui max_d = 0;
    for(auto e : add_vec){
        memset(dist, -1, sizeof(int)*n);
        dist[e] = 0;
        Q.push(e);
        
        while (!Q.empty()) {
            ui v = Q.front();
            Q.pop();
            if(dist[v] > max_d)
                max_d = dist[v];
            for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
                ui w = G0_edges[i];
                if(inVI[w] && dist[w] == -1){
                    dist[w] = dist[v] + 1;
                    Q.push(w);
                }
            }
        }//while
        
    }//for
    delete [] dist;
    
    if(max_d >= ubD)
        del_v_in_VI_ubD = true;
}*/

inline bool v1_in_v2Nei(ui v1, ui v2)
{
    int len = G0_x[v2];
    if(len == 0) return false;
    
    int s = pstart[v2];
    int e = pstart[v2] + len - 1;
    while (s <= e) {
        int mid = s + (e-s)/2;
        if(G0_edges[mid] == v1){
            return true;}
        if(G0_edges[mid] > v1){
            e = mid - 1;
        }
        else{
            s = mid + 1;
        }
    }
    return false;
}

void find_domS(vector<pair<ui, ui>> & domS)
{
    vector<pair<ui, ui>> vec_dv; //pair<degVIVR, vid>
    for(auto e : VIVR){
        if(inVR[e]){
            vec_dv.push_back(make_pair(degVIVR[e], e));
        }
    }
    if(vec_dv.size()==0){
        return;
    }
    sort(vec_dv.begin(), vec_dv.end(), less<>()); //increasing order
    
    bool * slt = new bool[n];
    memset(slt, 0, sizeof(bool)*n);
    
    for(auto e : vec_dv){
        ui v = e.second;
        if(slt[v]) continue;
        
        ui min_deg = INF;
        ui u = INF;
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            if(inVR[G0_edges[i]] || inVI[G0_edges[i]]){
                if(degVIVR[G0_edges[i]] < min_deg){
                    u = G0_edges[i];
                    min_deg = degVIVR[G0_edges[i]];
                }
            }
        }
        
        if(u == INF) continue;
        
        for(ui i = pstart[u]; i < pstart[u]+G0_x[u]; i++){
            ui w = G0_edges[i];
            if(inVR[w]){
                if(w == v || slt[w]) continue;
                bool flag = true;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui s = G0_edges[j];
                    if(inVR[s] || inVI[s]){
                        if(!v1_in_v2Nei(w,s) && w != s){
                            flag = false;
                            break;
                        }
                    }
                }
                if(flag){
                    slt[w] = 1;
                    slt[v] = 1;
                    domS.push_back(make_pair(w, v));
                    break;
                }
            }
        }
    }
    delete [] slt;
}

void find_domS_from_NEI(vector<pair<ui, ui>> & domS)
{
    if(NEI.size()==0 || NEI.size()==1) return;
    vector<pair<ui, ui>> vec_dv; //pair<degVIVR,vid>
    for(auto e : NEI)
        if(inVR[e]) vec_dv.push_back(make_pair(degVIVR[e], e));
    sort(vec_dv.begin(), vec_dv.end(), less<>());
    
    bool * slt = new bool[n];
    memset(slt, 0, sizeof(bool)*n);
    
    for(auto e : vec_dv){
        ui v = e.second;
        if(slt[v]) continue;
        
        ui min_deg = INF;
        ui u = INF;
        
        for(ui i = pstart[v]; i < pstart[v]+G0_x[v]; i++){
            if(inVR[G0_edges[i]] || inVI[G0_edges[i]]){
                if(degVIVR[G0_edges[i]] < min_deg){
                    u = G0_edges[i];
                    min_deg = degVIVR[G0_edges[i]];
                }
            }
        }
        
        if(u == INF) continue;
        
        for(ui i = pstart[u]; i < pstart[u]+G0_x[u]; i++){
            ui w = G0_edges[i];
            if(inNEI[w] > 0){
                if(w == v || slt[w]) continue;
    
                bool flag = true;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui s = G0_edges[j];
                    if(inVR[s] || inVI[s]){
                        if(!v1_in_v2Nei(w,s) && w != s){
                            flag = false;
                            break;
                        }
                    }
                }
                if(flag){
                    slt[w] = 1;
                    slt[v] = 1;
                    domS.push_back(make_pair(w, v));
                    break;
                }
            }
        }
        if(domS.size() > domS_Threshold)
            break;
    }
    delete [] slt;
}

void Cnm(int nn, int outLen, int starIdx, int mm, int * A, int AIdx)
{
    if(mm == 0){
        vector<ui> t_vec;
        for(int i = 0; i < outLen; i++){
//            cout<<A[i]<<",";
            t_vec.push_back(A[i]);
        }
        combs.push_back(t_vec);
//        cout<<endl;
        return;
    }
    int endIdx = nn-mm+1;
    for(int i = starIdx; i < endIdx; i++){
        A[AIdx] = i;
        Cnm(nn, outLen, i+1, mm-1, A, AIdx+1);
    }
    
}

void get_combs(ui dom_n)
{
    int * A = new int[dom_n];
    memset(A, 0, sizeof(int)*dom_n);
    for(ui i = 0; i <= dom_n; i++){
        Cnm(dom_n, i, 0, i, A, 0);
    }
    delete [] A;
}

void materialize_domS(const vector<pair<ui, ui>> & domS, vector<pair<vector<ui>, vector<ui>>> & M_domS)
{
    ui dom_n = (ui)domS.size();
    combs.clear();
    get_combs(dom_n);
    bool * in = new bool[dom_n];
    for(auto e : combs){
        vector<ui> t_parameter1;
        vector<ui> t_parameter2;
        memset(in, 0, sizeof(bool)*dom_n);
        for(auto x : e){
            in[x] = 1;
        }
        for(ui i = 0; i < dom_n; i++){
            if(in[i]){
                t_parameter1.push_back(domS[i].first);
                t_parameter1.push_back(domS[i].second);
//                t_parameter2.push_back(domS[i].first);
//                t_parameter2.push_back(domS[i].second);
            }
            else{
                t_parameter2.push_back(domS[i].second);
            }
        }
        M_domS.push_back(make_pair(t_parameter1, t_parameter2));
    }
    delete [] in;
}
    
void find_domS_of_ustar(ui ustar, vector<pair<double, ui>> & domS)
{
    for(auto e : NEI){
        if( inNEI[e] != 0 && degVI[e] <= degVI[ustar] && degVIVR[e] <= degVIVR[ustar] && e != ustar){
            ui ustar_sidx = pstart[ustar];
            ui ustar_eidx = pstart[ustar] + G0_x[ustar];
            bool be_dom = true;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i];
                if( (inVR[w] || inVI[w]) && w != ustar){
                    while(G0_edges[ustar_sidx] < w && ustar_sidx < ustar_eidx)
                        ++ ustar_sidx;
                    if(G0_edges[ustar_sidx] == w) continue;
                    else{
                        be_dom = false;
                        break;
                    }
                }
            }
            if(be_dom){
                double its_score = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]] && degVI[G0_edges[i]] != 0)
                      its_score += (double) 1/degVI[G0_edges[i]];
                }
                its_score += (double)degVIVR[e]/dMAX;
                domS.push_back(make_pair(its_score, e));
            }
        }
        if(domS.size() > domS_Threshold) break;
    }
    if(domS.size() > 1) sort(domS.begin(), domS.end(), greater<>());
}
    
void BB_dom_ustar()
{
    if(over_time_flag) return;
    
    double DurTime = (double)clock() / CLOCKS_PER_SEC - StartTime;
    
    if(DurTime > MaxTime)
        over_time_flag = true;
    
    // stop case for size constraint 
    if(VI.size() > h) return;

    // check to see if this instance min degree is bigger
    // if it is, update results and calculate new upD
    if(VI.size() >= l && VI.size() <= h){
        ui cur_min_deg = INF;
        for(auto e : VI){
            if(degVI[e] < cur_min_deg)
                cur_min_deg = degVI[e];
        }
        if(cur_min_deg > kl){
            kl = cur_min_deg;
            H = VI;
            for(ui d = 1; d <= h; d++){
                if(d == 1 || d == 2){
                    if(kl + d > h){
                        ubD = d-1;
                        break;
                    }
                }
                else{
                    ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                    if(h < min_n){
                        ubD = d - 1;
                        break;
                    }
                }
            }
        }
    }

 
    unordered_set<ui> new2VI;
    if(EXE_new2VI){
        for(auto e : VI){
            if(degVIVR[e] == kl+1){
                vector<ui> its_nei;
                for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w]){
                        its_nei.push_back(w);
                    }
                }
                if(its_nei.size() != 0){
                    for(auto x : its_nei){
                        new2VI.insert(x);
                    }
                }
            }
        }
        for(auto e : new2VI){
            if(inVR[e]){
                inVI[e] = 1;
                inVR[e] = 0;
                VI.push_back(e);
        //        cout<<e<<endl;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++ ){
                    if(inVI[G0_edges[i]]){
                        ++ degVI[G0_edges[i]];
                        ++ degVI[e];
                    }
                }
            }
        }
        if(VI.size() > h){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
            return;
        }
    } 
    
    // again check for new optimal // updatre ubD
    if(VI.size() >= l && VI.size() <= h){
        ui cur_min_deg = INF;
        for(auto e : VI){
            if(degVI[e] < cur_min_deg)
                cur_min_deg = degVI[e];
        }
        if(cur_min_deg > kl){
            kl = cur_min_deg;
            H = VI;
            for(ui d = 1; d <= h; d++){
                if(d == 1 || d == 2){
                    if(kl + d > h){
                        ubD = d - 1;
                        break;
                    }
                }
                else{
                    ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                    if(h < min_n){
                        ubD = d - 1;
                        break;
                    }
                }
            }
        }
    }
    if(VI.size() == h){
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
//                        cout<<"--degVI["<<w<<"] = "<<degVI[w]<<endl;
//                        cout<<"--degVI["<<e<<"] = "<<degVI[e]<<endl;
                    }
                }
            }
        }
        return;
    }
    
    NEI.clear();
    memset(inNEI, 0, sizeof(ui)*n);
    for(auto e : VI){
        for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
            if(inVR[G0_edges[i]]){
                if(inNEI[G0_edges[i]] == 0){
                    NEI.push_back(G0_edges[i]);
                    inNEI[G0_edges[i]] = 1;
                }
                else{
                    ++ inNEI[G0_edges[i]];
                }
            }
        }
    }

    vector<ui> del_from_VR;
    if(EXE_del_from_VR){
        for(auto e : NEI){
            if(inNEI[e] < kl+1){
                int lack = kl + 1 - inNEI[e];
                int bugt = h - (int)VI.size() - 1;
                if( lack > bugt ){
                    del_from_VR.push_back(e);
                    inVR[e] = 0;
                    inNEI[e] = 0;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            -- degVIVR[w];
                            -- degVIVR[e];
    //                                cout<<" -- degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
    //                                cout<<" -- degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                        }
                    }
                }
            }
        }
    }

    int ustar = -1;
    ustar = find_ustar_mindeg();
    
    if(ustar < 0){
        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
                    }
                }
            }
        }

        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
        }
        return;
    }

    vector<pair<double, ui>> domS;
    find_domS_of_ustar(ustar, domS);
    
    if(!domS.empty()){
        ++ domBr;
        for(ui round = 0; round < domS.size(); round++){
            ui domv = domS[round].second;
            
            VI.push_back(ustar);
            inVI[ustar] = 1;
            inVR[ustar] = 0;
            for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                if(inVI[G0_edges[i]]){
                    ++ degVI[ustar];
                    ++ degVI[G0_edges[i]];
                }
            }
            VI.push_back(domv);
            inVI[domv] = 1;
            inVR[domv] = 0;
            for(ui i = pstart[domv]; i < pstart[domv]+G0_x[domv]; i++){
                if(inVI[G0_edges[i]]){
                    ++ degVI[domv];
                    ++ degVI[G0_edges[i]];
                }
            }
            
            vector<ui> rmv_each_round;
            for(ui i = 0; i < round; i++){
                ui dv = domS[i].second;
                inVR[dv] = 0;
                rmv_each_round.push_back(dv);
                for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                    if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                        -- degVIVR[G0_edges[j]];
                        -- degVIVR[dv];
                    }
                }
            }
            
            vector<ui> rVR;
//            vector<bool> rVIVRbit;
            bool del_v_in_VI = false;
        //    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
            if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVR, rmv_each_round);
            if(EXE_core_maintenance){
                if(del_v_in_VI){
                    for(ui i = 0; i < rVR.size(); i++){
                        ui v = rVR[i];
                        inVR[v] = 1;
                        for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                            ui w = G0_edges[j];
                            if(inVI[w] || inVR[w]){
                                ++ degVIVR[w];
                                ++ degVIVR[v];
                            }
                        }
                    }
                    
                    for(ui i = 0; i < round; i++){
                        ui dv = domS[i].second;
                        inVR[dv] = 1;
                        for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                            if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                                ++ degVIVR[G0_edges[j]];
                                ++ degVIVR[dv];
                            }
                        }
                    }
                    
                    VI.pop_back();
                    inVI[domv] = 0;
                    inVR[domv] = 1;
                    for(ui i = pstart[domv]; i < pstart[domv]+G0_x[domv]; i++){
                        if(inVI[G0_edges[i]]){
                            -- degVI[domv];
                            -- degVI[G0_edges[i]];
                        }
                    }
                    
                    VI.pop_back();
                    inVI[ustar] = 0;
                    inVR[ustar] = 1;
                    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                        if(inVI[G0_edges[i]]){
                            -- degVI[ustar];
                            -- degVI[G0_edges[i]];
                        }
                    }
                    continue;
                }
            }
            
            ui ub = compute_ub();

            if(ub > kl){
                BB_dom_ustar();
            }

            if(EXE_core_maintenance){
                for(ui i = 0; i < rVR.size(); i++){
                    ui v = rVR[i];
                    inVR[v] = 1;
                    for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                        ui w = G0_edges[j];
                        if(inVI[w] || inVR[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[v];
                        }
                    }
                }
            }
            
            for(ui i = 0; i < round; i++){
                ui dv = domS[i].second;
                inVR[dv] = 1;
                for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                    if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                        ++ degVIVR[G0_edges[j]];
                        ++ degVIVR[dv];
                    }
                }
            }
            VI.pop_back();
            inVI[domv] = 0;
            inVR[domv] = 1;
            for(ui i = pstart[domv]; i < pstart[domv]+G0_x[domv]; i++){
                if(inVI[G0_edges[i]]){
                    -- degVI[domv];
                    -- degVI[G0_edges[i]];
                }
            }
            VI.pop_back();
            inVI[ustar] = 0;
            inVR[ustar] = 1;
            for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                if(inVI[G0_edges[i]]){
                    -- degVI[ustar];
                    -- degVI[G0_edges[i]];
                }
            }
        }
        
        for(ui i = 0; i < domS.size(); i++){
            ui dv = domS[i].second;
            inVR[dv] = 0;
            for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                    -- degVIVR[G0_edges[j]];
                    -- degVIVR[dv];
                }
            }
        }
        
        vector<ui> rVR;
//        vector<bool> rVIVRbit;
        bool del_v_in_VI = false;
    //    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
        if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVR, domS);
        if(EXE_core_maintenance){
            if(del_v_in_VI){
                for(ui i = 0; i < rVR.size(); i++){
                    ui v = rVR[i];
                    inVR[v] = 1;
                    for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                        ui w = G0_edges[j];
                        if(inVI[w] || inVR[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[v];
                        }
                    }
                }
                
                for(ui i = 0; i < domS.size(); i++){
                    ui dv = domS[i].second;
                    inVR[dv] = 1;
                    for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                        if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                            ++ degVIVR[G0_edges[j]];
                            ++ degVIVR[dv];
                        }
                    }
                }
                
                if(EXE_del_from_VR){
                    for(auto e : del_from_VR){
                        inVR[e] = 1;
                        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                            ui w = G0_edges[i];
                            if(inVR[w] || inVI[w]){
                                ++ degVIVR[w];
                                ++ degVIVR[e];
            //                    cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
            //                    cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                            }
                        }
                    }
                }
                if(EXE_new2VI){
                    for(auto e : new2VI){
                        VI.pop_back();
                        inVI[e] = 0;
                        inVR[e] = 1;
                        for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                            ui w = G0_edges[i];
                            if(inVI[w]){
                                -- degVI[w];
                                -- degVI[e];
                            }
                        }
                    }
                }
                return;
            }
        }
        
        ui ub = compute_ub();
        if(ub > kl){
            BB_dom_ustar();
        }
        
        if(EXE_core_maintenance){
            for(ui i = 0; i < rVR.size(); i++){
                ui v = rVR[i];
                inVR[v] = 1;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui w = G0_edges[j];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[v];
                    }
                }
            }
        }
        
        for(ui i = 0; i < domS.size(); i++){
            ui dv = domS[i].second;
            inVR[dv] = 1;
            for(ui j = pstart[dv]; j < pstart[dv]+G0_x[dv]; j++){
                if(inVI[G0_edges[j]] || inVR[G0_edges[j]]){
                    ++ degVIVR[G0_edges[j]];
                    ++ degVIVR[dv];
//                    cout<<"++ degVIVR["<<G0_edges[j]<<"] = "<<degVIVR[G0_edges[j]]<<endl;
//                    cout<<"++ degVIVR["<<dv<<"] = "<<degVIVR[dv]<<endl;
                }
            }
        }

        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
                    }
                }
            }
        }
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
        }
    }
    
    else{
        ++ binBr;
        VI.push_back(ustar);
        inVI[ustar] = 1;
        inVR[ustar] = 0;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            if(inVI[G0_edges[i]]){
                ++ degVI[G0_edges[i]];
                ++ degVI[ustar];
            }
        }
        
        ui ub = compute_ub();
        if(ub > kl){
            BB_dom_ustar();
        }
        
        
        VI.pop_back();
        inVI[ustar] = 0;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            if(inVI[G0_edges[i]]){
                -- degVI[G0_edges[i]];
                -- degVI[ustar];
    //            cout<<"-- degVI["<<G0_edges[i]<<"] = "<<degVI[G0_edges[i]]<<endl;
    //            cout<<"-- degVI["<<ustar<<"] = "<<degVI[ustar]<<endl;
            }
            if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
                -- degVIVR[G0_edges[i]];
                -- degVIVR[ustar];
    //            cout<<"-- degVIVR["<<G0_edges[i]<<"] = "<<degVIVR[G0_edges[i]]<<endl;
    //            cout<<"-- degVIVR["<<ustar<<"] = "<<degVIVR[ustar]<<endl;
            }
        }
        
        vector<ui> rVR;
//        vector<bool> rVIVRbit;
        bool del_v_in_VI = false;
    //    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
        if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVR, ustar);

        if(EXE_core_maintenance){
            if(del_v_in_VI){
                for(ui i = 0; i < rVR.size(); i++){
                    ui v = rVR[i];
                    inVR[v] = 1;
                    for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                        ui w = G0_edges[j];
                        if(inVI[w] || inVR[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[v];
                        }
                    }
                }
                inVR[ustar] = 1;
                for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[ustar];
                    }
                }
                
                for(auto e : del_from_VR){
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[e];
                        }
                    }
                }
                
                for(auto e : new2VI){
                    VI.pop_back();
                    inVI[e] = 0;
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVI[w]){
                            -- degVI[w];
                            -- degVI[e];
                        }
                    }
                }
                
                return;
            }
        }
        
        ub = compute_ub();
        if(ub > kl){
            BB_dom_ustar();
        }
        
        if(EXE_core_maintenance){
            for(ui i = 0; i < rVR.size(); i++){
                ui v = rVR[i];
                inVR[v] = 1;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui w = G0_edges[j];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[v];
                    }
                }
            }
        }
        
        inVR[ustar] = 1;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            ui w = G0_edges[i];
            if(inVI[w] || inVR[w]){
                ++ degVIVR[w];
                ++ degVIVR[ustar];
//                cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//                cout<<" ++ degVIVR["<<ustar<<"]="<<degVIVR[ustar]<<endl;
            }
        }
        
        
        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
                    }
                }
            }
        }
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
        }
    }
}


/*
defined & never used
void BB_dom_pairs()
{
    if(over_time_flag) return;
    
    double DurTime = (double)clock() / CLOCKS_PER_SEC - StartTime;
    
    if(DurTime > MaxTime)
        over_time_flag = true;
    
    if(VI.size() > h) return;
    if(VI.size() >= l && VI.size() <= h){
        ui cur_min_deg = INF;
        for(auto e : VI){
            if(degVI[e] < cur_min_deg)
                cur_min_deg = degVI[e];
        }
        if(cur_min_deg > kl){
            kl = cur_min_deg;
            H = VI;
            for(ui d = 1; d <= h; d++){
                if(d == 1 || d == 2){
                    if(kl + d > h){
                        ubD = d-1;
                        break;
                    }
                }
                else{
                    ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                    if(h < min_n){
                        ubD = d - 1;
                        break;
                    }
                }
            }
        }
    }
    if(VI.size() == h) return;

    vector<ui> new2VI;
    if(EXE_new2VI){
        for(auto e : VI){
            if(degVIVR[e] == kl+1){
                vector<ui> its_nei;
                for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w]){
                        its_nei.push_back(w);
                    }
                }
                if(degVI[e]+its_nei.size() == kl+1 && its_nei.size() != 0){
                    for(auto x : its_nei){
                        new2VI.push_back(x);
                    }
                }
            }
        }
        for(auto e : new2VI){
            inVI[e] = 1;
            inVR[e] = 0;
            VI.push_back(e);
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++ ){
                if(inVI[G0_edges[i]]){
                    ++ degVI[G0_edges[i]];
                    ++ degVI[e];
                }
            }
        }
    }
    
    NEI.clear();
    memset(inNEI, 0, sizeof(ui)*n);
    for(auto e : VI){
        for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
            if(inVR[G0_edges[i]]){
                if(inNEI[G0_edges[i]] == 0){
                    NEI.push_back(G0_edges[i]);
                    inNEI[G0_edges[i]] = 1;
                }
                else{
                    ++ inNEI[G0_edges[i]];
                }
            }
        }
    }
    
    vector<ui> del_from_VR;
    if(EXE_del_from_VR){
        for(auto e : NEI){
            if(inNEI[e] < kl+1){
                int lack = kl + 1 - inNEI[e];
                int bugt = h - (int)VI.size() - 1;
                if( lack > bugt ){
                    del_from_VR.push_back(e);
                    inVR[e] = 0;
                    inNEI[e] = 0;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            -- degVIVR[w];
                            -- degVIVR[e];
    //                                cout<<" -- degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
    //                                cout<<" -- degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                        }
                    }
                }
            }
        }
    }

    vector<pair<ui, ui>> domS;
    find_domS_from_NEI(domS);

    
    if(!domS.empty()){
        ++ domBr;
        vector<pair<vector<ui>, vector<ui>>> M_domS;
        materialize_domS(domS, M_domS);

        for(auto InOutPair : M_domS){
            vector<ui> add_to_VI = InOutPair.first;
            vector<ui> rmv_fr_VR = InOutPair.second;

            if(add_to_VI.size() > h - VI.size())
                continue;
            for(auto e : add_to_VI){
                VI.push_back(e);
                inVI[e] = 1;
                inVR[e] = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]]){
                        ++ degVI[G0_edges[i]];
                        ++ degVI[e];
//                        cout<<" ++ degVI["<<G0_edges[i]<<"]="<<degVI[G0_edges[i]]<<endl;
//                        cout<<" ++ degVI["<<e<<"]="<<degVI[e]<<endl;
                    }
                }
            }
            
            for(auto e : rmv_fr_VR){
                inVR[e] = 0;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
                        -- degVIVR[G0_edges[i]];
                        -- degVIVR[e];
//                        cout<<" -- degVIVR["<<G0_edges[i]<<"]="<<degVIVR[G0_edges[i]]<<endl;
//                        cout<<" -- degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                    }
                }
            }
            
                vector<ui> rVR;
//                vector<bool> rVIVRbit;
                bool del_v_in_VI = false;
            //    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
                if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVR, rmv_fr_VR);

                if(EXE_core_maintenance){
                    if(del_v_in_VI){
                        
                        for(ui i = 0; i < rVR.size(); i++){
                            ui v = rVR[i];
                            inVR[v] = 1;
                            for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                                ui w = G0_edges[j];
                                if(inVI[w] || inVR[w]){
                                    ++ degVIVR[w];
                                    ++ degVIVR[v];
                                }
                            }
                        }
                        
                        for(auto e : rmv_fr_VR){
                            inVR[e] = 1;
                            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                                if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
                                    ++ degVIVR[G0_edges[i]];
                                    ++ degVIVR[e];
                                }
                            }
                        }
                        
                        for(auto e : add_to_VI){
                            VI.pop_back();
                            inVI[e] = 0;
                            inVR[e] = 1;
                            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                                if(inVI[G0_edges[i]]){
                                    -- degVI[G0_edges[i]];
                                    -- degVI[e];
                                }
                            }
                        }
                        continue;
                    }
                }
            
            ui ub = compute_ub();
            if(ub > kl){
                BB_dom_pairs();
            }
            
            if(EXE_core_maintenance){
                for(ui i = 0; i < rVR.size(); i++){
                    ui v = rVR[i];
                    inVR[v] = 1;
                    for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                        ui w = G0_edges[j];
                        if(inVI[w] || inVR[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[v];
                        }
                    }
                }
            }
            
            for(auto e : rmv_fr_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
                        ++ degVIVR[G0_edges[i]];
                        ++ degVIVR[e];
//                        cout<<" ++ degVIVR["<<G0_edges[i]<<"]="<<degVIVR[G0_edges[i]]<<endl;
//                        cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                    }
                }
            }
    
            for(auto e : add_to_VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    if(inVI[G0_edges[i]]){
                        -- degVI[G0_edges[i]];
                        -- degVI[e];
//                        cout<<" -- degVI["<<G0_edges[i]<<"]="<<degVI[G0_edges[i]]<<endl;
//                        cout<<" -- degVI["<<e<<"]="<<degVI[e]<<endl;
                    }
                }
            }
        }
        
        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
    //                    cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
    //                    cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                    }
                }
            }
        }
        
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
        }
    }
    
    else{
        ++ binBr;
        int ustar = -1;
        if(srch_ord == 1) ustar = find_ustar();
        else if(srch_ord == 2) ustar = find_ustar_2phase();
        
        if(ustar < 0){
            if(EXE_del_from_VR){
                for(auto e : del_from_VR){
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[e];
    //                        cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
    //                        cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                        }
                    }
                }
            }
            
            if(EXE_new2VI){
                for(auto e : new2VI){
                    VI.pop_back();
                    inVI[e] = 0;
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVI[w]){
                            -- degVI[w];
                            -- degVI[e];
                        }
                    }
                }
            }
            
            return;
        }
        VI.push_back(ustar);
        inVI[ustar] = 1;
        inVR[ustar] = 0;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            if(inVI[G0_edges[i]]){
                ++ degVI[G0_edges[i]];
                ++ degVI[ustar];
            }
        }
        
        ui ub = compute_ub();
        if(ub > kl){
            BB_dom_pairs();
        }
        
        VI.pop_back();
        inVI[ustar] = 0;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            
            if(inVI[G0_edges[i]]){
                -- degVI[G0_edges[i]];
                -- degVI[ustar];
    //            cout<<"-- degVI["<<G0_edges[i]<<"] = "<<degVI[G0_edges[i]]<<endl;
    //            cout<<"-- degVI["<<ustar<<"] = "<<degVI[ustar]<<endl;
            }
            if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
                -- degVIVR[G0_edges[i]];
                -- degVIVR[ustar];
    //            cout<<"-- degVIVR["<<G0_edges[i]<<"] = "<<degVIVR[G0_edges[i]]<<endl;
    //            cout<<"-- degVIVR["<<ustar<<"] = "<<degVIVR[ustar]<<endl;
            }
        }
        
        vector<ui> rVIVR;
        vector<bool> rVIVRbit;
        bool del_v_in_VI = false;
    //    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
        if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVIVR, ustar);

        if(EXE_core_maintenance){
            if(del_v_in_VI){
                
                for(ui i = 0; i < rVIVR.size(); i++){
                    ui v = rVIVR[i];
                    if(rVIVRbit[i])
                        inVI[v] = 1;
                    else
                        inVR[v] = 1;
                    for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                        ui w = G0_edges[j];
                        if(inVI[w] || inVR[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[v];
                        }
                    }
                }
                
                inVR[ustar] = 1;
                for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[ustar];
                    }
                }
                
                for(auto e : del_from_VR){
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[e];
                        }
                    }
                }
                
                for(auto e : new2VI){
                    VI.pop_back();
                    inVI[e] = 0;
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVI[w]){
                            -- degVI[w];
                            -- degVI[e];
                        }
                    }
                }
                
                return;
            }
        }
        
        ub = compute_ub();
        if(ub > kl){
            BB_dom_pairs();
        }
        
        if(EXE_core_maintenance){
            for(ui i = 0; i < rVIVR.size(); i++){
                ui v = rVIVR[i];
                if(rVIVRbit[i])
                    inVI[v] = 1;
                else
                    inVR[v] = 1;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui w = G0_edges[j];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[v];
                    }
                }
            }
        }
        
        inVR[ustar] = 1;
        for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
            ui w = G0_edges[i];
            if(inVI[w] || inVR[w]){
                ++ degVIVR[w];
                ++ degVIVR[ustar];
            }
        }
        

        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
    //                    cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
    //                    cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                    }
                }
            }
        }
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
                    }
                }
            }
        }
    }
}*/

void BB()
{
    if(over_time_flag) return;
    
    double DurTime = (double)clock() / CLOCKS_PER_SEC - StartTime;
    
    if(DurTime > MaxTime) over_time_flag = true;
    

    // stop branching when size constraint broken
    // no need to check for feasability
    if(VI.size() > h) return; 
    
    // check against optimal
    if(VI.size() >= l && VI.size() <= h){ 
        ui cur_min_deg = INF;

        // find min degree of current partial solution to compare
        // against current best
        for(auto e : VI){
            if(degVI[e] < cur_min_deg)
                cur_min_deg = degVI[e];
        }

        // check if partial solution > current best
        if(cur_min_deg > kl){ 
            kl = cur_min_deg;
            H = VI;
            // search for up bound on diameter
            for(ui d = 1; d <= h; d++){
                if(d == 1 || d == 2){
                    if(kl + d > h){
                        ubD = d - 1;
                        break;
                    }
                }
                else{
                    ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                    if(h < min_n){
                        ubD = d - 1;
                        break;
                    }
                }
            }
        }
    }
    
    // VI doesn't shrink, it can't grow, stop checking
    if(VI.size() == h) return;
    
    /*
    for all critical vertices in partial solution
    add all adjacent vertices that are in candidate set
    to the partial solution(VI)

    if too many would be added (ie VI > h), do nothing
    */
    // RR2 flag, RR3 of paper
    unordered_set<ui> new2VI;
    if(EXE_new2VI){
        // for all vertices e in partial solution
        for(auto e : VI){
        // take only from promising
            if(degVIVR[e] == kl+1){
                
                vector<ui> its_nei;
                // add all neighbors of e 
                // (that are in candidate set) to new2VI
                for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w]){
                        its_nei.push_back(w);
                    }
                }

                if(its_nei.size() != 0){
                    for(auto x : its_nei){
                        new2VI.insert(x);
                    }
                }
            }
        }

        /* add all veritces in new2VI to VI*/

        // for all vertices e in new2VI
        for(auto e : new2VI){
            if(inVR[e]){
                // add e to VI, remove from VR
                inVI[e] = 1;
                inVR[e] = 0;
                VI.push_back(e);
                // update e and e's neighbor's degrees
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++ ){
                    if(inVI[G0_edges[i]]){
                        ++ degVI[G0_edges[i]];
                        ++ degVI[e];
                    }
                }
            }
        }


        /* if we added too many,
        revert back, and 
        remove all the ones we added */

        if(VI.size() > h){
            // for all e in new2VI (C U R)
            for(auto e : new2VI){
                // remove nodes in new2VI from VI
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;

                // fix degrees
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
//                        cout<<"--degVI["<<w<<"] = "<<degVI[w]<<endl;
//                        cout<<"--degVI["<<e<<"] = "<<degVI[e]<<endl;
                    }
                }
            }
            return;
        }
    }
    
    // check if (possibly...) new community is 
    // feasible & if better than current optimal
    if(VI.size() >= l && VI.size() <= h){
        ui cur_min_deg = INF;
        for(auto e : VI){
            if(degVI[e] < cur_min_deg)
                cur_min_deg = degVI[e];
        }
        if(cur_min_deg > kl){
            kl = cur_min_deg;
            H = VI;
            for(ui d = 1; d <= h; d++){
                if(d == 1 || d == 2){
                    if(kl + d > h){
                        ubD = d - 1;
                        break;
                    }
                }
                else{
                    ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                    if(h < min_n){
                        ubD = d - 1;
                        break;
                    }
                }
            }
        }
    }

    /*
    remove new vertices that were
    added from crticical vertices
    */
    // RR2 flag, RR3 of paper
    if(VI.size() == h){
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
//                        cout<<"--degVI["<<w<<"] = "<<degVI[w]<<endl;
//                        cout<<"--degVI["<<e<<"] = "<<degVI[e]<<endl;
                    }
                }
            }
        }
        return;
    }

    // degree of vertex v in candidate set in v U partial soluton
    NEI.clear();
    // degree 
    memset(inNEI, 0, sizeof(ui)*n);

    // RR3, RR1 FROM PAPER
    // populates NEI with vertices in candidate set that are
    // 1hop from a vertex in partial solution
    // inNEI is array holding deg(v) in partial solution U v
    for(auto e : VI){
        // for all neighbors v of e 
        for(ui i = pstart[e]; i < pstart[e] + G0_x[e]; i++){
            if(inVR[G0_edges[i]]){
                if(inNEI[G0_edges[i]] == 0){
                    NEI.push_back(G0_edges[i]);
                    inNEI[G0_edges[i]] = 1;
                }
                else{
                    ++ inNEI[G0_edges[i]];
                }
            }
        }
    }
    

    // reduction rule 1 of paper
    vector<ui> del_from_VR;
    if(EXE_del_from_VR){
        // check each vertex in candidate set
        for(auto e : NEI){
            // checking reduction rule
            // don't take min though, 
            // assume arbitrarily large degree (h - size:partial solution)
            if(inNEI[e] < kl+1){
                int lack = kl + 1 - inNEI[e];
                int bugt = h - (int)VI.size() - 1;
                if( lack > bugt ){
                    del_from_VR.push_back(e);
                    inVR[e] = 0;
                    inNEI[e] = 0;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            -- degVIVR[w];
                            -- degVIVR[e]; // removable line ? 
//                            cout<<" -- degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//                            cout<<" -- degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                        }
                    }
                }
            }
        }
    }

    int ustar = -1;
    if(srch_ord == 1) ustar = find_ustar();
    else if(srch_ord == 2) ustar = find_ustar_2phase();
    else if(srch_ord == 3) ustar = find_ustar_mindeg();
    else if(srch_ord == 4) ustar = find_ustar_link();
    else if(srch_ord == 5) ustar = find_ustar_random();
    
    // RR rules gave no u* to branch to, restore and end branch//return
    if(ustar < 0){
        if(EXE_del_from_VR){
            for(auto e : del_from_VR){
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVR[w] || inVI[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[e];
//                        cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//                        cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                    }
                }
            }
        }
        
        if(EXE_new2VI){
            for(auto e : new2VI){
                VI.pop_back();
                inVI[e] = 0;
                inVR[e] = 1;
                for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                    ui w = G0_edges[i];
                    if(inVI[w]){
                        -- degVI[w];
                        -- degVI[e];
//                        cout<<"--degVI["<<w<<"] = "<<degVI[w]<<endl;
//                        cout<<"--degVI["<<e<<"] = "<<degVI[e]<<endl;
                    }
                }
            }
        }
        return;
    }
    

    // done rules start branching (enumerating all possible subgraphs)
    // generate 2: one where the partial has v*
    // another where it does not
    // remove v* from candidate set in both

    // add ustar to partial solution
    // remove from candidate set
    VI.push_back(ustar);
    inVI[ustar] = 1;
    inVR[ustar] = 0;
    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
        if(inVI[G0_edges[i]]){
            ++ degVI[G0_edges[i]];
            ++ degVI[ustar];
//            cout<<"++ degVI["<<G0_edges[i]<<"] = "<<degVI[G0_edges[i]]<<endl;
//            cout<<"++ degVI["<<ustar<<"] = "<<degVI[ustar]<<endl;
        }
    }
    ui ub = compute_ub();
    
    // ensure it can't be pruned by its upper bound
    // ensure it is in reach from query node
    if(ub > kl && q_dist[ustar] <= ubD){
        BB();
    }
    
    // start generating second branch with no v*
    // remove v* from VI and branch
    VI.pop_back();
    inVI[ustar] = 0;
    
    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
        if(inVI[G0_edges[i]]){
            -- degVI[G0_edges[i]];
            -- degVI[ustar];
//            cout<<"-- degVI["<<G0_edges[i]<<"] = "<<degVI[G0_edges[i]]<<endl;
//            cout<<"-- degVI["<<ustar<<"] = "<<degVI[ustar]<<endl;
        }
        if(inVI[G0_edges[i]] || inVR[G0_edges[i]]){
            -- degVIVR[G0_edges[i]];
            -- degVIVR[ustar];
//            cout<<"-- degVIVR["<<G0_edges[i]<<"] = "<<degVIVR[G0_edges[i]]<<endl;
//            cout<<"-- degVIVR["<<ustar<<"] = "<<degVIVR[ustar]<<endl;
        }
    }

    vector<ui> rVR;
//    vector<bool> rVIVRbit;
    bool del_v_in_VI = false;
//    upd_based_kl(del_v_in_VI, rVIVR, rVIVRbit, ustar);
    
    if(EXE_core_maintenance) core_maintenance(del_v_in_VI, rVR, ustar);

    if(EXE_core_maintenance){
        if(del_v_in_VI){
            for(ui i = 0; i < rVR.size(); i++){
                ui v = rVR[i];
                inVR[v] = 1;
                for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                    ui w = G0_edges[j];
                    if(inVI[w] || inVR[w]){
                        ++ degVIVR[w];
                        ++ degVIVR[v];
//                        cout<<"++ degVIVR["<<w<<"] = "<<degVIVR[w]<<endl;
//                        cout<<"++ degVIVR["<<v<<"] = "<<degVIVR[v]<<endl;
                    }
                }
            }
            
            inVR[ustar] = 1;
            for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
                ui w = G0_edges[i];
                if(inVI[w] || inVR[w]){
                    ++ degVIVR[w];
                    ++ degVIVR[ustar];
                }
            }
            
            if(EXE_del_from_VR){
                for(auto e : del_from_VR){
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVR[w] || inVI[w]){
                            ++ degVIVR[w];
                            ++ degVIVR[e];
                        }
                    }
                }
            }
            
            if(EXE_new2VI){
                for(auto e : new2VI){
                    VI.pop_back();
                    inVI[e] = 0;
                    inVR[e] = 1;
                    for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                        ui w = G0_edges[i];
                        if(inVI[w]){
                            -- degVI[w];
                            -- degVI[e];
                        }
                    }
                }
            }
            return;
        }
    }
    
    ub = compute_ub();
        
    if(ub > kl){
        BB();
    }

    /*
         the rest of this method restores changes made to 
         VI the partial solution during reduction rules
    */


    if(EXE_core_maintenance){
        for(ui i = 0; i < rVR.size(); i++){
            ui v = rVR[i];
            inVR[v] = 1;
            for(ui j = pstart[v]; j < pstart[v]+G0_x[v]; j++){
                ui w = G0_edges[j];
                if(inVI[w] || inVR[w]){
                    ++ degVIVR[w];
                    ++ degVIVR[v];
//                    cout<<"++ degVIVR["<<w<<"] = "<<degVIVR[w]<<endl;
//                    cout<<"++ degVIVR["<<v<<"] = "<<degVIVR[v]<<endl;
                }
            }
        }
    }
    inVR[ustar] = 1;
    for(ui i = pstart[ustar]; i < pstart[ustar]+G0_x[ustar]; i++){
        ui w = G0_edges[i];
        if(inVI[w] || inVR[w]){
            ++ degVIVR[w];
            ++ degVIVR[ustar];
//            cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//            cout<<" ++ degVIVR["<<ustar<<"]="<<degVIVR[ustar]<<endl;
        }
    }
    
    if(EXE_del_from_VR){
        for(auto e : del_from_VR){
            inVR[e] = 1;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i];
                if(inVR[w] || inVI[w]){
                    ++ degVIVR[w];
                    ++ degVIVR[e];
//                    cout<<" ++ degVIVR["<<w<<"]="<<degVIVR[w]<<endl;
//                    cout<<" ++ degVIVR["<<e<<"]="<<degVIVR[e]<<endl;
                }
            }
        }
    }
    if(EXE_new2VI){
        for(auto e : new2VI){
            VI.pop_back();
            inVI[e] = 0;
            inVR[e] = 1;
            for(ui i = pstart[e]; i < pstart[e]+G0_x[e]; i++){
                ui w = G0_edges[i];
                if(inVI[w]){
                    -- degVI[w];
                    -- degVI[e];
//                    cout<<"-- degVI["<<w<<"] = "<<degVI[w]<<endl;
//                    cout<<"-- degVI["<<e<<"] = "<<degVI[e]<<endl;
                }
            }
        }
    }
}

void CSSC_BB()
{
    Timer timer;
    StartTime = (double)clock() / CLOCKS_PER_SEC;
    
    // compute core numbers, are in core[QID], 
    core_decomposition_linear_list();

    // UB of optimal
    ku = miv(core[QID], h-1);
    Timer t_for_heu;

    // compute heuristic 
    // for lower bound
    CSSC_heu();

    total_Heu_time += t_for_heu.elapsed();
    
    if(kl==ku){
        cout<<"heuristic find the OPT!"<<endl;
        cout<<"mindeg = "<<kl<<endl;
        cout<<"H.size = "<<H.size()<<endl;
        cout<<"time = "<<integer_to_string(timer.elapsed()).c_str()<<endl;
        return;
    }


    ubD = 0;
    // D is upper bounded by largest D s.t. n(kl, D) <= h
    // once we find a d that breaks this constraint, return
    // d - 1 for result
    // k = 1 is always h - 1, k=0 => bad graph input

    // all over the code this is posted, linear searching for D
    // can be improved via binary search since n(k, D) is a monotonic function
    if(kl<=1) ubD = h-1;
    else{
        for(ui d = 1; d <= h; d++){
            if(d == 1 || d == 2){
                if(kl + d > h){
                    ubD = d - 1;
                    break;
                }
            }
            else{
                ui min_n = kl + d + 1 + floor(d/3) * (kl - 2);
                if(h < min_n){
                    ubD = d - 1;
                    break;
                }
            }
        }
    }

    // if a vertex is farther from the query vertex than
    // the optimal upper bound diameter, we can
    // discard that vertex
    cal_query_dist();

    // remove cn(v) <= kl
    // resulting graph is in
    // p_start/G0/G0_x/G0_edges
    reduction_g();
        
    // VI : partial solution
    // VIVR: partial solution U candidate set
    VI.clear();
    VIVR.clear();
    
    // use to represent if a given vertex is 'in'
    // the partial solution or candidate set
    inVI = new bool[n];
    memset(inVI, 0, sizeof(bool)*n);
    inVR = new bool[n];
    memset(inVR, 0, sizeof(bool)*n);
    
    // degree of a vertex v in v U partial solution
    // degree of a vertex v in partial solution U candidate set
    degVI = new ui[n];
    memset(degVI, 0, sizeof(ui)*n);
    degVIVR = new ui[n];
    memset(degVIVR, 0, sizeof(ui)*n);

    // 
    inNEI = new ui[n];
    memset(inNEI, 0, sizeof(ui)*n);
    NEI_score = new double[n];
    memset(NEI_score, 0, sizeof(double)*n);
    

    // put all vertices in candidate set
    for(auto e : G0){
        VIVR.push_back(e); // (for all e in VR)
        inVR[e] = 1;        // for all v in G: if(inVR) 
        degVIVR[e] = G0_deg[e];
    }


    // start with required vertex in
    // partial solution 
    // remove it from candidate set
    VI.push_back(QID);
    inVI[QID] = 1;
    inVR[QID] = 0;
    
    over_time_flag = false;

    if(EXE_dom_ustar)
        BB_dom_ustar(); //with domination based branching
    else
        BB();           //baseline
    
    if(over_time_flag){
        cout<<"mindeg' = "<<kl<<endl;
        cout<<"H'.size = "<<H.size()<<endl;
        cout<<"time : "<<integer_to_string(timer.elapsed()).c_str()<<endl;
    }
    else{
        cout<<"mindeg = "<<kl<<endl;
        cout<<"H.size = "<<H.size()<<endl;
        cout<<"time = "<<integer_to_string(timer.elapsed()).c_str()<<endl;
    }
}

/*void graph_property()
defined never used again
{
    core_decomposition_linear_list();
    int max_core_number = 0;
    for(ui i = 0; i < n; i++){
        if(core[i] > max_core_number){
            max_core_number = core[i];
        }
    }
    
    int A[max_core_number+1];
    for(ui i=0; i <= max_core_number; i++){
        A[i] = 0;
    }
    for(ui i = 0; i<n; i++){
        ++ A[core[i]];
    }
    cout<<"core number distributiol: "<<endl;
    for(ui i=0; i <= max_core_number; i++){
        cout<<"core# = "<<i<<", total "<<A[i]<<" vertices."<<endl;
    }
    
    cout<<"core number distributioh: "<<endl;
    for(ui i=0; i <= max_core_number; i++){
        cout<<"core# >= "<<i;
        ui x = 0;
        for(ui j = i; j <= max_core_number; j++){
            x = x + A[j];
        }
        cout<<", total "<<x<<" vertices."<<endl;
    }
}*/

int main(int argc, const char * argv[]) {

    if(argc!=19){
        cout<<"wrong input parameters!"<<endl;exit(1);
    }
    
    l = atoi(argv[2]); //size LB
    h = atoi(argv[3]); //size UB
    QID = atoi(argv[4]); //Query vertex ID
    
    EXE_heu2 = atoi(argv[5]); //Heuristic strategy 1
    EXE_heu3 = atoi(argv[6]); //Heuristic strategy 2
    EXE_heu4 = atoi(argv[7]); //Heuristic strategy 3

    EXE_ub1 = atoi(argv[8]); //UB1
    EXE_ub2 = atoi(argv[9]); //UB2
    EXE_ub3 = atoi(argv[10]); //UB3
    EXE_ub3_optimization = atoi(argv[11]); //UB3 optimization
    
    EXE_core_maintenance = atoi(argv[12]); //reduction rule 1
    EXE_new2VI = atoi(argv[13]); //reduction rule 2
    EXE_del_from_VR = atoi(argv[14]); //reduction rule 3
    
    EXE_dom_ustar = atoi(argv[15]); //Dominating based branching rule
    domS_Threshold = atoi(argv[16]); //Dom pair threshold
    MaxTime = atoi(argv[17]); //OT
    srch_ord = atoi(argv[18]); //Branching order

    cout<<"Graph : "<<argv[1]<<", l = "<<l<<", h = "<<h<<", QID = "<<QID<<endl;
    cout<<"    Heu: "<<EXE_heu2<<","<<EXE_heu3<<","<<EXE_heu4;
    cout<<"    UBs: "<<EXE_ub1<<","<<EXE_ub2<<","<<EXE_ub3<<","<<EXE_ub3_optimization;
    cout<<"    Rdt: "<<EXE_core_maintenance<<","<<EXE_new2VI<<","<<EXE_del_from_VR;
    cout<<"    Dom: "<<EXE_dom_ustar<<","<<domS_Threshold;
    cout<<"    Tim: "<<MaxTime;
    cout<<"    Ord: "<<srch_ord;
    cout<<endl;

    load_graph(argv[1]);
    
    CSSC_BB();
    
    delete [] peel_sequence;
    delete [] degree;
    delete [] core;
    delete [] pstart;
    delete [] edges;
    delete [] q_dist;

    delete [] G0_edges;
    delete [] G0_x;
    delete [] G0_deg;
    
    delete [] inVI;
    delete [] inVR;
    delete [] degVI;
    
    delete [] inNEI;
    delete [] NEI_score;
    
    return 0;
}

