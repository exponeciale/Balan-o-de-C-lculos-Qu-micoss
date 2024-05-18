#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_ATOMOS 118
#define MAX_TERMOS 50

typedef struct {
    char termo[50];    // Termo da fórmula (por exemplo, "2H2O", "3CO2")
    char conjunto;     // Incógnita associada ao termo
} Associacao;

typedef struct {
    char symbol[3];  // Símbolo do átomo (por exemplo, "H", "O", "N", "Cl", "Al")
} Atomo;

void resetValores(int *numTermosEsquerda, int *numTermosDireita, int *j, int *ladoDireito);

void criarTabelaAssociacao(char *entrada, Associacao tabela[]) {
    const char *delimiters = "+=";
    char *token;
    int numTermos = 0;
    char incognita = 'a';
    // Tokenize a entrada usando os delimitadores
    token = strtok(entrada, delimiters);
    // Enquanto houver tokens
    while (token != NULL) {
        // Remove espaços em branco antes e depois do token
        char *termo = token;
        while (*termo == ' ') {
            termo++;
        }
        int length = strlen(termo);
        while (termo[length - 1] == ' ') {
            termo[length - 1] = '\0';
            length--;
        }
        // Salva o termo na tabela de associação e associa a incógnita correspondente
        strcpy(tabela[numTermos].termo, termo);
        tabela[numTermos].conjunto = incognita;
        // Próxima incógnita
        incognita++;
        // Próximo token
        token = strtok(NULL, delimiters);
        numTermos++;
    }
}

void extrairSimbolos(Associacao tabela[], int numTermos, Atomo atomos[], int *numAtomos) {
    *numAtomos = 0;
    for (int i = 0; i < numTermos; i++) {
        char *termo = tabela[i].termo;
        int j = 0;
        while (termo[j] != '\0') {
            // Verifica se o caractere é uma letra maiúscula
            if (termo[j] >= 'A' && termo[j] <= 'Z') {
                // Salva o símbolo de átomo como a letra maiúscula
                atomos[*numAtomos].symbol[0] = termo[j];
                // Se o próximo caractere é uma letra minúscula, salva também
                if (termo[j + 1] >= 'a' && termo[j + 1] <= 'z') {
                    atomos[*numAtomos].symbol[1] = termo[j + 1];
                    atomos[*numAtomos].symbol[2] = '\0'; // Terminador de string
                    (*numAtomos)++;
                    j++; // Avança mais um caractere
                } else {
                    // Se não houver letra minúscula, termina o símbolo de átomo
                    atomos[*numAtomos].symbol[1] = '\0';
                    (*numAtomos)++;
                }
                // Verifica se o símbolo já está na lista para evitar repetições
                for (int k = 0; k < *numAtomos - 1; k++) {
                    if (strcmp(atomos[*numAtomos - 1].symbol, atomos[k].symbol) == 0) {
                        (*numAtomos)--; // Remove o símbolo repetido
                    }
                }
            }
            j++;
        }
    }
}

void processaFormulaQuimica(char *formula, double coeficiente, double **tabela, Atomo *atomos, int numAtomos, int indiceTermo) {
    int len = strlen(formula);
    int i = 0;
    double coeficienteTermo = 1;  // Inicializa o coeficiente do termo como 1 por padrão
    // Verifica se o primeiro caractere da fórmula é um dígito
    if (isdigit(formula[i])) {
        // Usa sscanf para ler o coeficiente a partir da posição atual (i)
        sscanf(&formula[i], "%lf", &coeficienteTermo);
        // Atualiza o índice (i) para avançar após o coeficiente lido
        while (isdigit(formula[i]) && i < len) {
            i++;
        }
        // Multiplica coeficienteTermo pelo coeficiente total do termo
        coeficienteTermo *= coeficiente;
    }
    while (i < len) {
        if (isalpha(formula[i])) {
            // Início de um símbolo de átomo
            char simbolo[3] = { formula[i], '\0' };
            i++;
            while (islower(formula[i]) && i < len) {
                strncat(simbolo, &formula[i], 1);
                i++;
            }
            // Verifica se há um coeficiente numérico associado ao átomo
            double coeficienteAtomo = 1;
            if (isdigit(formula[i])) {
                sscanf(&formula[i], "%lf", &coeficienteAtomo);
                while (isdigit(formula[i]) && i < len) {
                    i++;
                }
            }
            // Encontra o índice do átomo na lista de átomos
            int indiceAtomo = -1;
            for (int k = 0; k < numAtomos; k++) {
                if (strcmp(atomos[k].symbol, simbolo) == 0) {
                    indiceAtomo = k;
                    break;
                }
            }
            // Preenche a tabela com a quantidade de átomos associada ao termo e desconhecido
            if (indiceAtomo != -1) {
                tabela[indiceAtomo][indiceTermo] += coeficiente * coeficienteTermo * coeficienteAtomo;
            }
        } else if (formula[i] == '(') {
            // Início de um grupo dentro de parênteses
            int j = i + 1;
            int profundidade = 1;
            // Encontra o fim do grupo dentro de parênteses
            while (j < len && profundidade > 0) {
                if (formula[j] == '(') {
                    profundidade++;
                } else if (formula[j] == ')') {
                    profundidade--;
                }
                j++;
            }
            // Processa o coeficiente do grupo dentro de parênteses, se houver
            double coeficienteGrupo = 1;
            if (j < len && isdigit(formula[j])) {
                sscanf(&formula[j], "%lf", &coeficienteGrupo);
                while (isdigit(formula[j]) && j < len) {
                    j++;
                }
            }
            // Processa recursivamente o grupo dentro de parênteses
            processaFormulaQuimica(&formula[i+1], coeficiente * coeficienteTermo * coeficienteGrupo, tabela, atomos, numAtomos, indiceTermo);
            i = j; // Atualiza o índice para após o grupo processado
        } else {
            i++;
        }
    }
}

double **imprimeTabelaIncognitas(Atomo *atomos, int numAtomos, Associacao *termos, int numTermos) {
    // Imprime cabeçalho com incognitas
    printf("\x1b[1m\x1b[33mTabela de Associação entre incognitas e Elementos:\x1b[0m\n");
    printf("\x1b[1m\x1b[33mX:\t\x1b[0m");
    for (int i = 0; i < numTermos; i++) {
        printf("\x1b[1m\x1b[33m%c\t\x1b[0m", termos[i].conjunto);
    }
    printf("\n");
    // Cria e inicializa a matriz para armazenar as quantidades associadas a cada incognitas
    double **tabela = malloc(numAtomos * sizeof(double*));
    for (int i = 0; i < numAtomos; i++) {
        tabela[i] = calloc(numTermos, sizeof(double));
    }
    // Preenche a tabela com as quantidades de cada átomo associadas a cada incognitas
    for (int j = 0; j < numTermos; j++) {
        processaFormulaQuimica(termos[j].termo, 1, tabela, atomos, numAtomos, j);
    }
    // Imprime os elementos associados a cada átomo para cada incognitas
    for (int i = 0; i < numAtomos; i++) {
        printf("\x1b[1m\x1b[33m%s:\t\x1b[0m", atomos[i].symbol);
        for (int j = 0; j < numTermos; j++) {
            printf("%.0lf\x1b[33m%c\t\x1b[0m", tabela[i][j], termos[j].conjunto);
        }
        printf("\n");
    }
    return tabela;
}

double **imprimirResultados(double **tabela, int numAtomos, int numTermos, Associacao *termos, char incognitaGlobal) {
    // Aloca memória para a matriz de doubles
    double **resultadosFormatados = (double **)malloc((numAtomos + 1) * sizeof(double *));
    for (int i = 0; i < numAtomos + 1; i++) {
        resultadosFormatados[i] = (double *)malloc((numTermos + 1) * sizeof(double));
    }
    // Preenche a primeira linha com o valor 1 na coluna associada à incógnita e 0 nas demais colunas
    for (int j = 0; j < numTermos; j++) {
        resultadosFormatados[0][j] = (termos[j].conjunto == incognitaGlobal) ? 1 : 0;
    }
    resultadosFormatados[0][numTermos] = 1; // Adiciona 1 para a nova coluna
    // Preenche as demais linhas com os valores da tabela e 0 para a nova coluna
    for (int i = 0; i < numAtomos; i++) {
        for (int j = 0; j < numTermos; j++) {
            resultadosFormatados[i + 1][j] = tabela[i][j];
        }
        resultadosFormatados[i + 1][numTermos] = 0; // Adiciona 0 para a nova coluna
    }
    // Retorna a matriz de doubles com os resultados formatados
    return resultadosFormatados;
}

void obterColunaMaiorSoma(double **tabela, int numAtomos, int numTermos, char *incognitaGlobal) {
    double maiorSoma = 0;
    char incognita = 'a';
    
    for (int j = 0; j < numTermos; j++) {
        double soma = 0;
        for (int i = 0; i < numAtomos; i++) {
            soma += tabela[i][j];
        }
        if (soma > maiorSoma) {
            maiorSoma = soma;
            *incognitaGlobal = incognita;
        }
        incognita++;
    }
}

void imprimirResultadosFormatados(double **resultadosFormatados, int numLinhas, int numColunas) {
    // Imprime os resultados formatados
    for (int i = 0; i < numLinhas; i++) {
        for (int j = 0; j < numColunas; j++) {
            printf("%.2lf ", resultadosFormatados[i][j]);
        }
        printf("\n");
    }
}

void printMatrix(double **matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j <= cols; j++) {
            printf("%.2f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
void rowOperations(double **matrix, int rows, int cols, int row1, int row2, double factor) {
    // Aplica a operação: row1 = row1 + factor * row2
    for (int j = 0; j <= cols; j++) {
        matrix[row1][j] += factor * matrix[row2][j];
    }
}

void makeIdentity(double **matrix, int rows, int cols) {
    int lead = 0;
    
    for (int r = 0; r < rows; r++) {
        if (lead >= cols) break;
        // Encontra uma linha com um valor não zero na coluna lead
        int i = r;
        while (fabs(matrix[i][lead]) < 1e-12) {
            i++;
            if (i == rows) {
                i = r;
                lead++;
                if (lead == cols) {
                    break;
                }
            }
        }
        // Troca a linha atual com a linha encontrada
        for (int j = 0; j <= cols; j++) {
            double temp = matrix[r][j];
            matrix[r][j] = matrix[i][j];
            matrix[i][j] = temp;
        }
        // Normaliza a linha atual
        double divisor = matrix[r][lead];
        if (fabs(divisor) > 1e-12) {
            for (int j = 0; j <= cols; j++) {
                matrix[r][j] /= divisor;
            }
        }
        // Zera todas as outras entradas na coluna de chumbo
        for (i = 0; i < rows; i++) {
            if (i != r) {
                double factor = -matrix[i][lead];
                rowOperations(matrix, rows, cols, i, r, factor);
            }
        }
        lead++;
    }
}
void printPenultimateColumn(double column[], int N) {
    printf("Penúltima coluna:\n");
    for (int i = 0; i < N; i++) {
        printf("%.2f\n", fabs(column[i]));
    }
}
int main() {
    char entrada[100];
    Associacao tabela[100]; // Array de estruturas Associacao
    Atomo atomos[100];       // Array de estruturas Atomo
    int numAtomos;
    char incognitaGlobal;
    int numTermosEsquerda = 0;
    int numTermosDireita = 0;
    char termoAtual[20]; // Para armazenar temporariamente cada termo
    int j = 0;           // Índice para o termo atual
    int ladoDireito = 0; // Flag para indicar o lado direito da equação
    char continuar[4];   // Variável para armazenar a entrada do usuário

    do {
        printf("\x1b[1m\x1b[32mInforme a equação química:\x1b[0m\n");
        fgets(entrada, sizeof(entrada), stdin);
        entrada[strcspn(entrada, "\n")] = '\0';
        int comprimento = strlen(entrada);
        int i = 0;
        int numTermos = 1; // Começa com 1 termo (o primeiro elemento antes do sinal de igualdade)
        char entrada_original[100];
        char entrada_original2[100];
        strcpy(entrada_original, entrada);
        strcpy(entrada_original2, entrada);

        while (i < comprimento) {
            if (entrada[i] == '+' || entrada[i] == '=') {
                // Encontrou um sinal de adição ou igualdade, indicando um novo termo na equação
                numTermos++;
            }
            // Avança para o próximo caractere
            i++;
        }
        // Exibir a quantidade de termos na equação
        printf("\n");
        // Cria a tabela de associação
        criarTabelaAssociacao(entrada, tabela);
        // Extrai os símbolos dos átomos da equação
        extrairSimbolos(tabela, numTermos, atomos, &numAtomos);
        // Chamada da função imprimeTabelaIncognitas
        double **resultadoTabela = imprimeTabelaIncognitas(atomos, numAtomos, tabela, numTermos);
        obterColunaMaiorSoma(resultadoTabela, numAtomos, numTermos, &incognitaGlobal);
        // Uso do resultadoTabela para operações posteriores
        double **resultadoFormatacao = imprimirResultados(resultadoTabela, numAtomos, numTermos, tabela, incognitaGlobal);
        makeIdentity(resultadoFormatacao, numAtomos + 1, numTermos + 1);
        printf("\x1b[1m\x1b[31m\nMatriz identidade:\x1b[0m\n");
        printMatrix(resultadoFormatacao, numAtomos + 1, numTermos + 1);
        // Verifica se a última linha é nula
        int isNull = 1;
        for (int j = 0; j < numTermos + 1; j++) {
            if (fabs(resultadoFormatacao[numAtomos][j]) > 1e-12) {
                isNull = 0;
                break;
            }
        }
        if (isNull) {
            // Isola a penúltima coluna das linhas
            double penultimateColumn[numAtomos];
            for (int i = 0; i < numAtomos; i++) {
                penultimateColumn[i] = resultadoFormatacao[i][numTermos];
            }
            int comprimento2 = strlen(entrada_original);
            while (i < comprimento2) {
                if (entrada_original[i] == '+' || entrada_original[i] == '=') {
                    // Encontrou um sinal de adição ou igualdade, indicando um novo termo na equação
                    numTermos++;
                }
                // Avança para o próximo caractere
                i++;
            }
            // Reescreve a equação e conta os termos novamente usando a entrada original
            // Separando os termos antes e depois da igualdade
            printf("\x1b[1m\x1b[36mReescrevendo a equação balanceada:\x1b[0m\n");
            for (int i = 0; i < comprimento2; i++) {
                if (entrada_original[i] == '+' || entrada_original[i] == '=') {
                    termoAtual[j] = '\0'; // Terminando a string atual
                    if (strlen(termoAtual) > 0) {
                        if (ladoDireito) {
                            printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s", fabs(penultimateColumn[numTermosDireita+numTermosEsquerda]), termoAtual);
                            numTermosDireita++;
                            // Adicionando sinal de '+' se não for o último termo do lado    direito
                            if (numTermosDireita < comprimento2 && entrada_original[i] != '\0')
                                printf(" \x1b[1m\x1b[36m+\x1b[0m ");
                        } else {
                            printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s", fabs(penultimateColumn[numTermosEsquerda]), termoAtual);
                            numTermosEsquerda++;
                            // Adicionando sinal de '+' se não for o último termo do lado    esquerdo
                            if (numTermosEsquerda < comprimento2 && entrada_original[i] != '=')
                            printf(" \x1b[1m\x1b[36m+\x1b[0m ");
                            else if (numTermosEsquerda < comprimento2)
                            printf(" \x1b[1m\x1b[36m=\x1b[0m ");
                        }
                        j = 0; // Reiniciando o índice do termo atual
                    }
                    if (entrada_original[i] == '=') {
                        ladoDireito = 1; // Indicando que estamos no lado direito da equação
                    }
                } else if (entrada_original[i] != ' ') {
                    termoAtual[j] = entrada_original[i];
                    j++;
                }
            }
            // Imprimir o último termo do lado direito, se houver
            if (ladoDireito) {
                termoAtual[j] = '\0'; // Terminando a string atual
                if (strlen(termoAtual) > 0) {
                    printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s", fabs(penultimateColumn[numTermosDireita+numTermosEsquerda]), termoAtual);
                }
            }
        } else {
            // Isolar a última coluna de todas as linhas até que sejam todos inteiros
            int multiplier = 1;
            int allInteger = 0;
            while (!allInteger) {
                allInteger = 1;
                for (int i = 0; i < numAtomos; i++) {
                    double result = resultadoFormatacao[i][numTermos] * multiplier;
                    if (fabs(result - round(result)) > 1e-12) {
                        allInteger = 0;
                        break;
                    }
                }
                if (!allInteger) {
                    multiplier++;
                }
            }
            // Multiplica a última coluna pelo multiplicador encontrado
            for (int i = 0; i <= numAtomos; i++) {
                resultadoFormatacao[i][numTermos] *= multiplier;
            }    
            // Armazenar os valores modificados da última coluna
            int comprimento3 = strlen(entrada_original2);
            while (i < comprimento3) {
                if (entrada_original2[i] == '+' || entrada_original2[i] == '=') {
                    // Encontrou um sinal de adição ou igualdade, indicando um novo termo na equação
                    numTermos++;
                }
                // Avança para o próximo caractere
                i++;
            }
            double modifiedLastColumn[numAtomos + 1];
            for (int i = 0; i <= numAtomos; i++) {
            	modifiedLastColumn[i] = fabs(resultadoFormatacao[i][numTermos]);
            }
            // Reescrever a equação original adicionando os valores modificados da última coluna
            printf("\x1b[1m\x1b[36mReescrevendo a equação balanceada:\x1b[0m\n");
            for (int i = 0; i < comprimento3; i++) {
            	if (entrada_original2[i] == '+' || entrada_original2[i] == '=') {
        	        termoAtual[j] = '\0'; // Terminando a string atual
   	                 if (strlen(termoAtual) > 0) {
               	        if (ladoDireito) {
                            printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s",modifiedLastColumn[numTermosDireita+numTermosEsquerda], termoAtual);
            	            // Adicionando o valor modificado da última coluna ao termo atual
            	            printf(" \x1b[1m\x1b[36m+\x1b[0m ");
                            numTermosDireita++;       
                	    } else {
                            printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s", modifiedLastColumn[numTermosEsquerda], termoAtual);
                            numTermosEsquerda++;
                            // Adicionando sinal de '+' se não for o último termo do lado esquerdo
                            if (numTermosEsquerda < comprimento3 && entrada_original2[i] != '=')
                            printf(" \x1b[1m\x1b[36m+\x1b[0m ");
                            else if (numTermosEsquerda < comprimento3)
                            printf(" \x1b[1m\x1b[36m+\x1b[0m ");
                	    }
                        j = 0; // Reiniciando o índice do termo atual
                    }
                    if (entrada_original2[i] == '=') {
                        ladoDireito = 1; // Indicando que estamos no lado direito da equação
                    }
                } else if (entrada_original2[i] != ' ') {
                    termoAtual[j] = entrada_original2[i];
                    j++;
                }
            }
            // Imprimir o último termo do lado direito, se houver
            if (ladoDireito) {
                termoAtual[j] = '\0'; // Terminando a string atual
                if (strlen(termoAtual) > 0) {
                    printf("\x1b[1m\x1b[36m%.0f\x1b[0m %s\n", termoAtual, modifiedLastColumn[numTermosDireita+numTermosEsquerda]);
                }
            }
        }
        // Liberar a memória alocada
        for (int i = 0; i < numAtomos; i++) {
            free(resultadoTabela[i]);
        }
        free(resultadoTabela);
	/*
        // Liberar a memória alocada para cada linha da matriz
        for (int i = 0; i < numAtomos + 1; i++) {
            free(resultadoFormatacao[i]);
        }
	    // Liberar a memória alocada para a matriz
	    free(resultadoFormatacao);
        */
        printf("\n\n\x1b[1m\x1b[33mDeseja continuar (\x1b[0m\x1b[1m\x1b[32msim\x1b[0m/\x1b[1m\x1b[31mfim\x1b[0m\x1b[1m\x1b[33m)?\x1b[0m ");
        fgets(continuar, sizeof(continuar), stdin);
        getchar(); // Limpa o buffer do teclado

        if (strcmp(continuar, "sim\n") == 0 || strcmp(continuar, "sim") == 0) {
        // Limpar a tela
        #ifdef _WIN32
            system("cls");
        #else
            system("clear");
        #endif
        resetValores(&numTermosEsquerda, &numTermosDireita, &j, &ladoDireito);
    }
} 
while (strcmp(continuar, "fim\n") != 0 && strcmp(continuar, "fim") != 0);
return 0;
}
void resetValores(int *numTermosEsquerda, int *numTermosDireita, int *j, int *ladoDireito) {
    *numTermosEsquerda = 0;
    *numTermosDireita = 0;
    *j = 0;
    *ladoDireito = 0;
}
