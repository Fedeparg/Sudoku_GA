#include <ga/GASimpleGA.h> //  Algoritmo Genetico simple
#include <ga/GA1DArrayGenome.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Definimos la estructura plantilla que mantendrá los valores fijos.
struct plantilla
{
    int tam;
    int *fijo;
};

float Objective(GAGenome &); // Funcion objetivo --> al final
GABoolean Termina(GAGeneticAlgorithm &); // Funcion de terminacion --> al final
void leerSudoku(struct plantilla *S, char *nombreF);
void imprimirResultado(GASimpleGA ga, int tam);
void InicioSudoku(GAGenome& g);
int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2);
int MutacionSudoku(GAGenome& g,float pmut);

int main(int argc, char **argv)
{

// Declaramos variables para los parametros del GA y las inicializamos

    char * sudoku_file (argv[1]);
    int popsize = atoi(argv[2]);
    int selector = atoi(argv[3]);
    float pcross = atof(argv[4]);
    float pmut = atof(argv[5]);
    int ngen = 12000;

    plantilla plantilla;


    leerSudoku(&plantilla, sudoku_file);

// Conjunto enumerado de alelos --> valores posibles de cada gen del genoma

    GAAlleleSet<int> alelos;
    for(int i=0; i<plantilla.tam; i++)
        alelos.add(i);

// Creamos el genoma y definimos operadores de inicio, cruce y mutaci�n

    GA1DArrayAlleleGenome<int> genome(plantilla.tam*plantilla.tam,alelos,Objective,&plantilla);
    genome.initializer(InicioSudoku);
    genome.crossover(CruceSudoku);
    genome.mutator(MutacionSudoku);

// Creamos el algoritmo genetico

    GASimpleGA ga(genome);

// Inicializamos - minimizar funcion objetivo, tama�o poblacion, n� generaciones,
// pr. cruce y pr. mutacion, selecci�n y le indicamos que evolucione.

    ga.minimaxi(-1);
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pCrossover(pcross);
    ga.pMutation(pmut);


    // Elegimos el selector según el parámetro indicado
    string selec = "";
    if (selector == 0){
        GATournamentSelector selector;
        ga.selector(selector);
        selec = "GATournament";
    } else if (selector == 1){
        GARouletteWheelSelector selector;
        ga.selector(selector);
        selec = "GARouletteWheel";
    }

    cout << "Fichero de entrada: " << sudoku_file << endl << endl;
    cout << "Parametros:    - Tamano poblacion: " << popsize << endl;
    cout << "               - Operador de seleccion: " << selec << endl;
    cout << "               - Numero de generaciones: " << ngen << endl;
    cout << "               - Probabilidad cruce: " << pcross << endl;
    cout << "               - Probabilidad mutacion: " << pmut << endl << endl;

    ga.terminator(Termina);
    ga.evolve(1);

    // Imprimimos el mejor individuo que encuentra el GA y su valor fitness
    imprimirResultado(ga, plantilla.tam);
}


// Funcion para inicializar el sudoku
void InicioSudoku(GAGenome& g)
{

    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;
    struct plantilla * plantilla1;
    plantilla1 = (struct plantilla *) genome.userData();

    int aux[plantilla1->tam]; // Generamos un array del tamaño del lado del sudoku


    for(int f=0; f<plantilla1->tam; f++)
    {

        for(int j=0; j<plantilla1->tam; j++)
            aux[j]=0; // Limpiamos el array con 0

        // Rellenamos de forma aleatoria el array
        for(int j=1; j<=plantilla1->tam; j++)
        {
            int v=GARandomInt(0,plantilla1->tam-1);
            while (aux[v]!=0)
                v=(v+1)%plantilla1->tam;
            aux[v]=j;
        }

        /**
         * A continuación comprobamos que los valores del array auxuliar no
         * entren en conflicto con los fijos del sudoku. Estos son reconocidos
         * por ser distintos de 0
        **/
        int i=0;

        while(i<plantilla1->tam)
        {

            while((plantilla1->fijo[(f*plantilla1->tam)+i]==0) && (i<plantilla1->tam))
                i++;

            if (i<plantilla1->tam)
            {

                bool encontrado=false;
                for(int j=0; (j<plantilla1->tam) && (!encontrado); j++)
                    if (aux[j]==plantilla1->fijo[(f*plantilla1->tam)+i])
                    {
                        encontrado=true;
                        aux[j]=aux[i];
                    }

                aux[i]=plantilla1->fijo[(f*plantilla1->tam)+i];
            }
            i++;

        }

        // Volcamos el array auxiliar al cromosoma
        for(int c=0; c<plantilla1->tam; c++)
            genome.gene((f*plantilla1->tam)+c,aux[c]);
    }
}

// Funcion para leer el sudoku de un fichero
void leerSudoku(struct plantilla *S,char *nombreF)
{
    // Abrimos el fichero
    ifstream f(nombreF);

    // Almacenamos el tamano
    f>>S->tam;

    /* Creamos un array del tamano del sudoku que contendra
    ** los valores fijos del sudoku */
    S->fijo = new int[S->tam*S->tam];

    // Llenamos el array con los valores
    for(int i=0; i<S->tam*S->tam; i++)
        f>>S->fijo[i];

    // Cerramos el fichero
    f.close();
}



bool checkColumna(int col[], int * check, int tam)
{
    bool repe=false;

    for(int i=0; i<tam; i++)
        check[i]=0;

    for(int i=0; i<tam; i++)
        check[col[i]-1]++;
    for(int i=0; i<tam; i++)
        if (check[i]>1)
            repe=true;

    return repe;
}


// Función de mutación
int MutacionSudoku(GAGenome& g,float pmut)
{

    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;

    struct plantilla * plantilla1;
    plantilla1 = (struct plantilla *) genome.userData();
    int nmut=0;
    int aux;
    int fil;
    bool fila;

    int caux[plantilla1->tam];
    int *checkC=new int[plantilla1->tam];

    if (pmut<=0.0)
        return 0;

    /** Primero se comprueba si el valor que se está analizando es parte del 
     * sudoku original en caso de no serlo, se modificará o no el valor según 
     * la probabilidad de mutación
    **/
    for(int f=0; f<plantilla1->tam; f++)
        for(int c=0; c<plantilla1->tam; c++)
            if (plantilla1->fijo[(f*plantilla1->tam)+c]==0)
            {
                if (GAFlipCoin(pmut) )
                {
                    // Si se decide que el valor debe mutar, se determina si
                    // mutar filas o columnas con probabilidad 0.5
                    if (GAFlipCoin(0.5))
                        fila = true;
                    else
                        fila = false;

                    /**
                     * En el caso de ser por columnas, se comprueba si hay
                     * repetidos en la columna actual y si el resultado es
                     * negativo, no se ejecuta mutación alguna.
                    **/
                    if (!fila)
                    {

                        for(int j=0; j<plantilla1->tam; j++)
                            caux[j]=genome.gene((j*plantilla1->tam)+c);

                        /**
                         * Si hay repetidos, se selecciona un número con más de 
                         * una repetición y otro que no aparezca
                        **/
                        if (checkColumna(caux,checkC,plantilla1->tam))
                        {
                            int v1 = GARandomInt(0,plantilla1->tam-1);
                            while (checkC[v1]<=1)
                                v1=(v1+1)%plantilla1->tam;
                            v1++;
                            int v2 = GARandomInt(0,plantilla1->tam-1);
                            while (checkC[v2]!=0)
                                v2=(v2+1)%plantilla1->tam;
                            v2++;
                            /**
                             * Como desconocemos las posiciones de los elementos
                             * repetidos, debemos buscarlas.
                            **/
                            bool encontrado = false;
                            for(int j=0; j<plantilla1->tam && !encontrado; j++)
                                if ((plantilla1->fijo[j*(plantilla1->tam)+c]==0)&&(genome.gene(j*(plantilla1->tam)+c)==v1))
                                {
                                    // Cuando la encontremos, cambiaos su valor
                                    // por la del número sin apariciones
                                    encontrado = true;
                                    genome.gene((j*plantilla1->tam)+c,v2);
                                    fil = j;
                                }

                            int col=(c+1)%plantilla1->tam;
                            while(genome.gene((fil*plantilla1->tam)+col)!=v2)
                                col=(col+1)%plantilla1->tam;
                            if (plantilla1->fijo[(fil*plantilla1->tam)+col]==0)
                            {
                                nmut++;
                                genome.gene((fil*plantilla1->tam)+col,v1);
                            }
                            else
                            {
                                genome.gene((fil*plantilla1->tam)+c,v1);
                            }

                        }

                    }
                    // EN el caso de ser por filas, las posiciones de dos
                    // valores pertenecientes a la misma fila, son permutados
                    else
                    {
                        int v1 = (c + 1) %plantilla1->tam;
                        while ((plantilla1->fijo[(f*plantilla1->tam)+v1]!=0))
                            v1=(v1+1)%plantilla1->tam;
                        aux = genome.gene((f*plantilla1->tam)+c);
                        genome.gene((f*plantilla1->tam)+c,genome.gene((f*plantilla1->tam)+v1));
                        genome.gene((f*plantilla1->tam)+v1,aux);
                        nmut++;
                    }
                }
            }

    return nmut;
}


// Funcion de cruce 
int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2)
{
    /**
     * Es una operación de cruce por 1 punto que genera dos hijos
    **/

    GA1DArrayAlleleGenome<int> & m = (GA1DArrayAlleleGenome<int> &)p1;
    GA1DArrayAlleleGenome<int> & p = (GA1DArrayAlleleGenome<int> &)p2;


    struct plantilla * plantilla1 = (struct plantilla *) m.userData();
    int n=0;

    // Seleccionamos un punto aleatorio
    int punto1=GARandomInt(0,m.length());
    while ((punto1%plantilla1->tam)!=0)
        punto1++;
        // Otro punto aleatorio calculado en base al primero
    int punto2=m.length()-punto1;

    /**
     * Generamos dos hijos copiando hasta el primer punto de un padre y el resto
     * del otro (representado por el valor dle punto 2)
    **/
    if (c1)
    {
        GA1DArrayGenome<int> & h1 = (GA1DArrayGenome<int> &) *c1;
        h1.copy(m,0,0,punto1); // el metodo copy esta definido en la clase GA1DArrayGenome
        h1.copy(p,punto1,punto1,punto2);
        n++;
    }

    if (c2)
    {

        GA1DArrayGenome<int> & h2 = (GA1DArrayGenome<int> &) *c2;
        h2.copy(p,0,0,punto1);
        h2.copy(m,punto1,punto1,punto2);
        n++;
    }

    return n;

}

// Fucnion que calcula el fitness de cada individuo
float Objective(GAGenome& g)
{
    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;
    float puntuacion = 0;
    int tamano = sqrt(genome.length());
    int subcuadrados = sqrt(tamano);

    // Recorremos tantas veces como números tenga el sudoku
    for (int numero = 1; numero<=tamano; numero++)
    {
        // Comprobamos por filas
        for(int i = 0; i <genome.length(); i++)
        {
            int contador = 0;
            if (genome.gene(i) == numero)
            {
                contador++;
            }

            // Si llegamos al final de una fila
            if (i%tamano == tamano-1)
            {
                // Si hemos contado mas de un numero similar
                if (contador > 1)
                {
                    // Aumentar puntuación
                    puntuacion += contador-1;
                }
                // Reseteamos contador para la siguiente fila
                contador = 0;
            }
        }

        // Recorremos por columnas
        for (int i = 0; i <tamano; i++)
        {
            int contador = 0;
            for (int j = i; j<genome.length(); j=j+tamano)
            {
                if (genome.gene(j) == numero)
                    contador++;
            }
            if (contador > 1)
                puntuacion += contador-1;
        }

        // Comprobamos los cuadrados
        for(int i = 0; i < genome.length(); i+=subcuadrados)
        {
            int contador = 0;
            int j = i;
            for(int k = 1; k <= tamano; k++)
            {
                if(genome.gene(j) == numero)
                {
                    contador++;
                }

                if (j%(tamano/subcuadrados) == subcuadrados-1)
                {
                    j+= tamano-subcuadrados;
                }
                j+=1;
            }

            if (contador > 1)
            {
                puntuacion+= contador-1;
            }
            contador = 0;

            if (i%tamano == tamano-subcuadrados)
            {
                i += (subcuadrados*subcuadrados*(subcuadrados-1));
            }
        }
    }
    return puntuacion;
}

// Condicion de terminacion
GABoolean Termina(GAGeneticAlgorithm & ga)
{
    // Si el fitness del minimo es 0 o hemos llegado al límite de generaciones
    if ( (ga.statistics().minEver()==0) ||
            (ga.statistics().generation()==ga.nGenerations()) )
        return gaTrue;
    else
        return gaFalse;
}


// Funcion para imprimir el resultado
void imprimirResultado(GASimpleGA ga, int tam)
{
    cout << "La solucion:" << endl;
    GA1DArrayAlleleGenome<int> & mejor = (GA1DArrayAlleleGenome<int> &)ga.statistics().bestIndividual();

    for(int i = 0; i < tam*tam; i++)
    {
        cout << mejor.gene(i) << " ";

        if(i%tam == tam-1)
        {
            cout << endl;
        }
    }
    cout << endl;

    cout << "Con valor fitness " << ga.statistics().minEver() << endl;
}
