#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>

using namespace std;

typedef unsigned long long ull; 
typedef long long ll;

unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine gen(seed);

ull random(ull low, ull high){
    uniform_int_distribution<ull> dist(low, high);
    return dist(gen);
}
int randomInteger(int low, int high){
    uniform_int_distribution<int> dist(low, high);
    return dist(gen);
}
double randomDouble(double low, double high){
    uniform_real_distribution<double> dist(low, high);
    return dist(gen);
}

void imprimeBits(ull valor, int numeroBits){
    for(int i = (numeroBits - 1); i >= 0; i--){
        ull aux = 1ULL << i;
        if(valor & aux)
            cout<<"1";
        else 
            cout<<"0";
        
    }
    cout << "\n";
}    

double F6(double x, double y){
    double resultado = 0.5 - (pow(sin(sqrt((x*x) +(y*y))), 2.0) - 0.5)/(pow(1.0 + 0.01 * ((x*x) +(y*y)), 2.0));   
    return resultado;
}

void geraFilhos(ull mae, ull pai, ull &filho1, ull &filho2, ull bitMaskLower, ull bitMaskUpper){
    filho1 = mae & bitMaskUpper;
    filho1 |= pai & bitMaskLower;
    
    filho2 = pai & bitMaskUpper;
    filho2 |= mae & bitMaskLower;
}

ull controleMutacao(ull individuo, vector<ull> &mascaras, double taxaMutacao){
    for(int i = 0; i < 44; i++){
        ull aleatorio = randomDouble(0.0, 1.0001);  
        if(aleatorio <= taxaMutacao){
            ull um = 1ULL << i;
            if(individuo & um)
                individuo &= mascaras[i];
            else
                individuo |= um;
        }
    }
    return individuo;
}
int melhorIndividuo(vector<ull> &populacao, ull bitMaskLower, double fator){
    ull x =  populacao[0] >> 22;
    ull y =  populacao[0] & bitMaskLower;  
    double xD = x * fator - 100.0; 
    double yD = y * fator - 100.0;
    double resultado = F6(xD, yD);
    double melhorSolucao = resultado;
    int indiceMelhorSolucao = 0;
    for(int i = 1; i < (int)populacao.size(); i++){
             x =  populacao[i] >> 22;
             y =  populacao[i] & bitMaskLower;  
             xD = x * fator - 100.0; 
             yD = y * fator - 100.0;
             resultado = F6(xD, yD);
            if(resultado > melhorSolucao){
                melhorSolucao = resultado;
                indiceMelhorSolucao = i;
            }
    }  
    return indiceMelhorSolucao;
    
}
int piorIndividuo(vector<ull> &populacao, ull bitMaskLower, double fator){
    ull x =  populacao[0] >> 22;
    ull y =  populacao[0] & bitMaskLower;  
    double xD = x * fator - 100.0; 
    double yD = y * fator - 100.0;
    double resultado = F6(xD, yD);
    double piorSolucao = resultado;
    int indicePiorSolucao = 0;
    for(int i = 1; i < (int)populacao.size(); i++){
             x =  populacao[i] >> 22;
             y =  populacao[i] & bitMaskLower;  
             xD = x * fator - 100.0; 
             yD = y * fator - 100.0;
             resultado = F6(xD, yD);
            if(resultado < piorSolucao){
                piorSolucao = resultado;
                indicePiorSolucao = i;
            }
         
    }
    return piorSolucao;
    
}
int main(){
    ull maxValue = 4194303ULL; // 2^22 - 1
    double fator = (double)(200.0/maxValue);
    double taxaCrossover =  0.65;
    double taxaMutacao = 0.008;
    int pop_size = 10000;
    int num_ger = 100;
    ull minValue = 0ULL;
    ull bitMaskLower = 0ULL;
    ull bitMaskUpper = 0ULL;
    ull bitMaskAux = 0ULL;
    vector<ull> mascarasMutacao; 
    
    for(int i = 44; i < 64; i++){
        bitMaskAux |= 1ULL << i;
    }
    bitMaskAux =  ~bitMaskAux;
    
    for(int i = 0; i < 44; i++){
        ull auxMascara = 0;
        auxMascara |= 1ULL << i;
        auxMascara = ~auxMascara; 
        auxMascara &= bitMaskAux;
        mascarasMutacao.push_back(auxMascara);
    }
    
    for(int i = 0; i < 22; i++){
        bitMaskLower |= 1ULL << i;
    }
    for(int i = 22; i < 44; i++){
        bitMaskUpper |= 1ULL << i;
    }
    
    //Gerar População Aleatória
    vector<ull> populacao(pop_size);
    for(int i = 0; i < pop_size; i++){
        ull numero = random(minValue, maxValue);
        ull gene = numero << 22;
    
        numero = random(minValue, maxValue);
        gene  |= numero;
        
        populacao[i] = gene;
    }   
    
    //Treino
    for(int i = 0; i < num_ger; i++){
        //Classificação dos indivíduos
        vector<double> acumulado(pop_size);
        for(int j = 0; j < pop_size; j++){
            ull x =  populacao[j] >> 22;
            ull y =  populacao[j] & bitMaskLower;  
            double xD = x * fator - 100.0; 
            double yD = y * fator - 100.0;
            double resultado = F6(xD, yD);
            if(j > 0) acumulado[j] = resultado + acumulado[j - 1];
            else acumulado[j] = resultado;
        }
        double maior = acumulado[pop_size - 1] + 0.1;
        //cout<<"maior acumulado: "<<acumulado[pop_size-1] << "\n";
        //cout<<"maior: " << maior << "\n";
        
        //Gerar a nova população
        vector<ull> novaPopulacao;
        for(int j = 0; j < (pop_size/2); j++){
            //Seleciona Pai
            double aleatoria1 = randomDouble(0.0, maior);
            //cout<<aleatoria1<<"\n";
            bool flag = false;
            int indicePai;
            for(int k = 0; k < pop_size; k++){
                if(acumulado[k] >= aleatoria1){
                    flag = true; 
                    indicePai = k;
                    break;
                }
            }
            if(!flag)
                indicePai = pop_size - 1;
            
            //Seleciona Mãe
            flag = false;
            double aleatoria2 = randomDouble(0.0, maior);
            //cout<<aleatoria2<<"\n";
            int indiceMae;
            for(int k = 0; k < pop_size; k++){
                if(acumulado[k] >= aleatoria2){
                    flag = true; 
                    indiceMae = k;
                    break;
                }
            }
            if(!flag)
                indiceMae = pop_size - 1;
            
            //Gera os filhos e os coloca na nova população
            double crossover = randomDouble(0.0, 1.0001);
            if(crossover <= taxaCrossover){
                    ull filho1 = 0ULL, filho2 = 0ULL;
                    geraFilhos(populacao[indiceMae], populacao[indicePai], filho1, filho2, bitMaskLower, bitMaskUpper);
                    novaPopulacao.push_back(filho1);
                    novaPopulacao.push_back(filho2);
            }else{
                    novaPopulacao.push_back(populacao[indicePai]);
                    novaPopulacao.push_back(populacao[indiceMae]);
            }   
            
        }  
        //Mutação
        for(int j = 0; j < pop_size;  j++){
           novaPopulacao[j] = controleMutacao(novaPopulacao[j], mascarasMutacao, taxaMutacao);
        }
        int indiceMelhorPai = melhorIndividuo(populacao, bitMaskLower, fator); 
        int indicePiorFilho = piorIndividuo(novaPopulacao, bitMaskLower, fator); 
        novaPopulacao[indicePiorFilho] = populacao[indiceMelhorPai];
        populacao = novaPopulacao;  
    }
    
    

    int indiceMelhorSolucao = melhorIndividuo(populacao, bitMaskLower, fator);
    ull x =  populacao[indiceMelhorSolucao] >> 22;
    ull y =  populacao[indiceMelhorSolucao] & bitMaskLower;  
    double xD = x * fator - 100.0; 
    double yD = y * fator - 100.0;
    double melhorSolucao = F6(xD, yD);
    
   
    cout <<"Melhor solução encontrada: " << melhorSolucao << "\n";
    cout << "Índice da melhor solução : "  << indiceMelhorSolucao << "\n";
    return 0;  
}    
