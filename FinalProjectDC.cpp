#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <omp.h>
#include <cmath>
#define N 16
using namespace std;

string codon[] = {
            "ata", "atc", "att", "atg",
            "aca", "acc", "acg", "act",
            "aac", "aat", "aaa", "aag",
            "agc", "agt", "aga", "agg",
            "cta", "ctc", "ctg", "ctt",
            "cca", "ccc", "ccg", "cct",
            "cac", "cat", "caa", "cag",
            "cga", "cgc", "cgg", "cgt",
            "gta", "gtc", "gtg", "gtt",
            "gca", "gcc", "gcg", "gct",
            "gac", "gat", "gaa", "gag",
            "gga", "ggc", "ggg", "ggt",
            "tca", "tcc", "tcg", "tct",
            "ttc", "ttt", "tta", "ttg",
            "tac", "tat", "taa", "tag",
            "tgc", "tgt", "tga", "tgg" };

char trancodon[] = { 'I','I','I','M','T','T','T','T','N','N',
            'K','K','S','S','R','R','L','L','L','L','P','P','P',
            'P','H','H','Q','Q','R','R','R','R','V','V','V','V','A','A','A','A','D','D',
            'E','E','G','G','G','G','S','S','S','S','F','F','L','L','Y','Y','_','_','C','C','_','W' };

string dna = "";

void read_seq()
{
    fstream newfile;

    char tpp[3000];
    newfile.open("InputSeq.txt", ios::in);
    if (newfile.is_open())
    {
        string tp = " ";
        while (getline(newfile, tp))
        {

            tp.erase(remove(tp.begin(), tp.end(), ' '), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '1'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '2'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '3'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '4'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '5'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '6'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '7'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '8'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '9'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '0'), tp.end());
            tp.erase(remove(tp.begin(), tp.end(), '/'), tp.end());

            dna.append(tp);
        }
        newfile.close();
    }
    cout << dna << "\n";

}

void serialcode()
{
    string protein = "", codo = "", s = "", amino = "";
    char I = 'I', M = 'M', T = 'T', n = 'N', K = 'K', S = 'S', R = 'R', L = 'L', P = 'P', H = 'H', Q = 'Q', V = 'V', A = 'A', D = 'D', E = 'E', G = 'G', F = 'F', Y = 'Y', stop = '_', C = 'C', W = 'W';
    char letters[21] = { 'I', 'M','T','N', 'K','S', 'R', 'L', 'P', 'H', 'Q', 'V', 'A', 'D','E','G', 'F','Y', '_','C', 'W' };
    int count = 0, ex = 0, en = 0, id = 0;
    int i, j, x, k, lengthP, doneboo = 0, index = 0;
    float percent = 0.0, sum = 0.0;
    double start, end;
    en = ((dna.size()) - (dna.size()) % 3) - 1;

    omp_set_num_threads(N);

    start = omp_get_wtime();

    for (int i = 0; i <= en; i += 3)
    {
        for (int j = i; j < (i + 3); j++)
        {

            s = dna[j];
            codo.append(s);
        }
        for (int x = 0; x < 64; x++)
        {
            if (codo == codon[x])
            {
                amino = trancodon[x];
                protein.append(amino);
            }
        }
        codo = "";
    }

    cout << protein << "\n";
    int counter[21];
    lengthP = protein.length();

    for (k = 0; k < 21; k++)
    {
        counter[k] = 0;
    }

    for (int i = 0; i < lengthP; i++)
    {
        if (protein[i] == I)
        {
            counter[0] = counter[0] + 1;
        }
        if (protein[i] == M)
        {
            counter[1] = counter[1] + 1;
        }
        if (protein[i] == T)
        {
            counter[2] = counter[2] + 1;
        }
        if (protein[i] == n)
        {
            counter[3] = counter[3] + 1;
        }
        if (protein[i] == K)
        {
            counter[4] = counter[4] + 1;
        }
        if (protein[i] == S)
        {
            counter[5] = counter[5] + 1;
        }
        if (protein[i] == R)
        {
            counter[6] = counter[6] + 1;
        }
        if (protein[i] == L)
        {
            counter[7] = counter[7] + 1;
        }
        if (protein[i] == P)
        {
            counter[8] = counter[8] + 1;
        }
        if (protein[i] == H)
        {
            counter[9] = counter[9] + 1;
        }
        if (protein[i] == Q)
        {
            counter[10] = counter[10] + 1;
        }
        if (protein[i] == V)
        {
            counter[11] = counter[11] + 1;
        }
        if (protein[i] == A)
        {
            counter[12] = counter[12] + 1;
        }
        if (protein[i] == D)
        {
            counter[13] = counter[13] + 1;
        }
        if (protein[i] == E)
        {
            counter[14] = counter[14] + 1;
        }
        if (protein[i] == G)
        {
            counter[15] = counter[15] + 1;
        }
        if (protein[i] == F)
        {
            counter[16] = counter[16] + 1;
        }
        if (protein[i] == Y)
        {
            counter[17] = counter[17] + 1;
        }
        if (protein[i] == stop)
        {
            counter[18] = counter[18] + 1;
        }
        if (protein[i] == C)
        {
            counter[19] = counter[19] + 1;
        }
        if (protein[i] == W)
        {
            counter[20] = counter[20] + 1;
        }
    }
    end = omp_get_wtime();
    for (int i = 0; i < 21; i++)
    {
        percent = (float)(((float)counter[i] / (float)lengthP) * 100);
        sum = sum + percent;
        cout << "Frequency of:  " << letters[i] << " is " << counter[i] << " Percentage is  " << percent << "%" << "\n";
    }

    cout << "Time" << end - start << "\n";
    cout << "Length of protein: " << lengthP << "\n";
    cout << "Total Frequency: " << sum << "%" << "\n";
}

void paralelcode()
{
    string protein = "", s = "", codo = "", amino = "";
    char I = 'I', M = 'M', T = 'T', n = 'N', K = 'K', S = 'S', R = 'R', L = 'L', P = 'P', H = 'H', Q = 'Q', V = 'V', A = 'A', D = 'D', E = 'E', G = 'G', F = 'F', Y = 'Y', stop = '_', C = 'C', W = 'W';
    char letters[21] = { 'I', 'M','T','N', 'K','S', 'R', 'L', 'P', 'H', 'Q', 'V', 'A', 'D','E','G', 'F','Y', '_','C', 'W' };
    int count = 0, ex = 0, en = 0, id = 0;
    int i, j, x, k, lengthP, doneboo = 0, index = 0;
    float percent = 0.0, sum = 0.0;
    double start, end;
    en = ((dna.size()) - (dna.size()) % 3) - 1;

    omp_set_num_threads(N);

    start = omp_get_wtime();
#pragma omp parallel for schedule(dynamic) firstprivate(s,dna,codo,amino) shared(en,protein)
    for (int i = 0; i <= en; i += 3)
    {
#pragma omp parallel for
        for (int j = i; j < (i + 3); j++)
        {

            s = dna[j];
            codo.append(s);
        }
#pragma omp parallel for
        for (int x = 0; x < 64; x++)
        {
            if (codo == codon[x])
            {
                amino = trancodon[x];
#pragma omp critical
                protein.append(amino);
            }
        }
        codo = "";
    }

    cout << protein << "\n";
    int counter[21];
    lengthP = protein.length();
#pragma omp parallel for shared(counter,k) 
    for (k = 0; k < 21; k++)
    {
        counter[k] = 0;
    }
#pragma omp parallel for
    for (int i = 0; i < lengthP; i++)
    {

        if (protein[i] == I)
        {
#pragma omp critical
            counter[0] = counter[0] + 1;
        }
        if (protein[i] == M)
        {
#pragma omp critical
            counter[1] = counter[1] + 1;
        }
        if (protein[i] == T)
        {
#pragma omp critical
            counter[2] = counter[2] + 1;
        }
        if (protein[i] == n)
        {
#pragma omp critical
            counter[3] = counter[3] + 1;
        }
        if (protein[i] == K)
        {
#pragma omp critical
            counter[4] = counter[4] + 1;
        }
        if (protein[i] == S)
        {
#pragma omp critical
            counter[5] = counter[5] + 1;
        }
        if (protein[i] == R)
        {
#pragma omp critical
            counter[6] = counter[6] + 1;
        }
        if (protein[i] == L)
        {
#pragma omp critical
            counter[7] = counter[7] + 1;
        }
        if (protein[i] == P)
        {
#pragma omp critical
            counter[8] = counter[8] + 1;
        }
        if (protein[i] == H)
        {
#pragma omp critical
            counter[9] = counter[9] + 1;
        }
        if (protein[i] == Q)
        {
#pragma omp critical
            counter[10] = counter[10] + 1;
        }
        if (protein[i] == V)
        {
#pragma omp critical
            counter[11] = counter[11] + 1;
        }
        if (protein[i] == A)
        {
#pragma omp critical
            counter[12] = counter[12] + 1;
        }
        if (protein[i] == D)
        {
#pragma omp critical
            counter[13] = counter[13] + 1;
        }
        if (protein[i] == E)
        {
#pragma omp critical
            counter[14] = counter[14] + 1;
        }
        if (protein[i] == G)
        {
#pragma omp critical
            counter[15] = counter[15] + 1;
        }
        if (protein[i] == F)
        {
#pragma omp critical
            counter[16] = counter[16] + 1;
        }
        if (protein[i] == Y)
        {
#pragma omp critical
            counter[17] = counter[17] + 1;
        }
        if (protein[i] == stop)
        {
#pragma omp critical
            counter[18] = counter[18] + 1;
        }
        if (protein[i] == C)
        {
#pragma omp critical
            counter[19] = counter[19] + 1;
        }
        if (protein[i] == W)
        {
#pragma omp critical
            counter[20] = counter[20] + 1;
        }
    }
    end = omp_get_wtime();
    for (int i = 0; i < 21; i++)
    {
        percent = (float)(((float)counter[i] / (float)lengthP) * 100);
        sum = sum + percent;
        cout << "Frequency of:  " << letters[i] << " is " << counter[i] << " Percentage is  " << percent << "%" << "\n";

    }

    cout << "Time" << end - start << "\n";
    cout << "Length of protein: " << lengthP << "\n";
    cout << "Total Frequency: " << sum << "%" << "\n";
}

void prallelSectionCode()
{

    string protein1 = "", protein2 = "", protein3 = "", protein4 = "", protein5 = "", protein6 = "", protein7 = "", protein8 = "", protein9 = "", protein10 = "", protein11 = "", protein12 = "", protein13 = "", protein14 = "", protein15 = "", protein16 = "";
    string codo = "", s = "", amino = "";
    char I = 'I', M = 'M', T = 'T', n = 'N', K = 'K', S = 'S', R = 'R', L = 'L', P = 'P', H = 'H', Q = 'Q', V = 'V', A = 'A', D = 'D', E = 'E', G = 'G', F = 'F', Y = 'Y', stop = '_', C = 'C', W = 'W';
    char letters[21] = { 'I', 'M','T','N', 'K','S', 'R', 'L', 'P', 'H', 'Q', 'V', 'A', 'D','E','G', 'F','Y', '_','C', 'W' };
    int count = 0, ex = 0, en = 0, id = 0;
    int i, j, x, k, lengthP, doneboo = 0, index = 0;
    float percent = 0.0, sum = 0.0;
    double start, end;
    en = ((dna.size()) - (dna.size()) % 3) - 1;

    omp_set_num_threads(N);

    start = omp_get_wtime();
#pragma omp parallel firstprivate(codo,amino,s,en)
    {
#pragma omp sections 
        {
#pragma omp section
            {
                for (int i = 0; i < floor((en / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein1.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor((en / 16)) + 1; i < floor(((2 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein2.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((2 * en) / 16)) + 1; i < floor(((3 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein3.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((3 * en) / 16)) + 2; i < floor(((4 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein4.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((4 * en) / 16)) + 2; i < floor(((5 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein5.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((5 * en) / 16)); i < floor(((6 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein6.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((6 * en) / 16)); i < floor(((7 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein7.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((7 * en) / 16)) + 1; i < floor(((8 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein8.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((8 * en) / 16)) + 1; i < floor(((9 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein9.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((9 * en) / 16)) + 1; i < floor(((10 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein10.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((10 * en) / 16)) + 2; i < floor(((11 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein11.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((11 * en) / 16)) + 2; i < floor(((12 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein12.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((12 * en) / 16)); i < floor(((13 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein13.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((13 * en) / 16)); i < floor(((14 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein14.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((14 * en) / 16)) + 1; i < floor(((15 * en) / 16)); i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein15.append(amino);
                        }
                    }
                    codo = "";
                }
            }
#pragma omp section
            {
                for (int i = floor(((15 * en) / 16)) + 1; i < (16 * en) / 16; i += 3)
                {
                    for (int j = i; j < (i + 3); j++)
                    {

                        s = dna[j];
                        codo.append(s);
                    }

                    for (int x = 0; x < 64; x++)
                    {
                        if (codo == codon[x])
                        {
                            amino = trancodon[x];
                            protein16.append(amino);
                        }
                    }
                    codo = "";
                }
            }
        }
    }
    protein1.append(protein2);
    protein1.append(protein3);
    protein1.append(protein4);
    protein1.append(protein5);
    protein1.append(protein6);
    protein1.append(protein7);
    protein1.append(protein8);
    protein1.append(protein9);
    protein1.append(protein10);
    protein1.append(protein11);
    protein1.append(protein12);
    protein1.append(protein13);
    protein1.append(protein14);
    protein1.append(protein15);
    protein1.append(protein16);

    cout << protein1 << "\n";
    int counter[21];
    lengthP = protein1.length();
    for (k = 0; k < 21; k++)
    {
        counter[k] = 0;
    }

#pragma omp parallel 
    {
#pragma omp sections
        {
#pragma omp section
            {
                for (int i = 0; i < floor(lengthP / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }

#pragma omp section
            {
                for (int i = floor(lengthP / 16); i < floor((2 * lengthP) / 16); i++)
                {
                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = floor((2 * lengthP) / 16); i < floor((3 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((3 * lengthP) / 16); i < floor((4 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((4 * lengthP) / 16); i < floor((5 * lengthP) / 16); i++)
                {
                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((5 * lengthP) / 16); i < floor((6 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((6 * lengthP) / 16); i < floor((7 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }

            }
#pragma omp section
            {
                for (int i = ((7 * lengthP) / 16); i < floor((8 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((8 * lengthP) / 16); i < floor((9 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((9 * lengthP) / 16); i < floor((10 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((10 * lengthP) / 16); i < floor((11 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((11 * lengthP) / 16); i < floor((12 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((12 * lengthP) / 16); i < floor((13 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((13 * lengthP) / 16); i < floor((14 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((14 * lengthP) / 16); i < floor((15 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
#pragma omp section
            {
                for (int i = ((15 * lengthP) / 16); i < floor((16 * lengthP) / 16); i++)
                {

                    if (protein1[i] == I)
                    {
#pragma omp critical
                        counter[0] = counter[0] + 1;
                    }
                    if (protein1[i] == M)
                    {
#pragma omp critical
                        counter[1] = counter[1] + 1;
                    }
                    if (protein1[i] == T)
                    {
#pragma omp critical
                        counter[2] = counter[2] + 1;
                    }
                    if (protein1[i] == n)
                    {
#pragma omp critical
                        counter[3] = counter[3] + 1;
                    }
                    if (protein1[i] == K)
                    {
#pragma omp critical
                        counter[4] = counter[4] + 1;
                    }
                    if (protein1[i] == S)
                    {
#pragma omp critical
                        counter[5] = counter[5] + 1;
                    }
                    if (protein1[i] == R)
                    {
#pragma omp critical
                        counter[6] = counter[6] + 1;
                    }
                    if (protein1[i] == L)
                    {
#pragma omp critical
                        counter[7] = counter[7] + 1;
                    }
                    if (protein1[i] == P)
                    {
#pragma omp critical
                        counter[8] = counter[8] + 1;
                    }
                    if (protein1[i] == H)
                    {
#pragma omp critical
                        counter[9] = counter[9] + 1;
                    }
                    if (protein1[i] == Q)
                    {
#pragma omp critical
                        counter[10] = counter[10] + 1;
                    }
                    if (protein1[i] == V)
                    {
#pragma omp critical
                        counter[11] = counter[11] + 1;
                    }
                    if (protein1[i] == A)
                    {
#pragma omp critical
                        counter[12] = counter[12] + 1;
                    }
                    if (protein1[i] == D)
                    {
#pragma omp critical
                        counter[13] = counter[13] + 1;
                    }
                    if (protein1[i] == E)
                    {
#pragma omp critical
                        counter[14] = counter[14] + 1;
                    }
                    if (protein1[i] == G)
                    {
#pragma omp critical
                        counter[15] = counter[15] + 1;
                    }
                    if (protein1[i] == F)
                    {
#pragma omp critical
                        counter[16] = counter[16] + 1;
                    }
                    if (protein1[i] == Y)
                    {
#pragma omp critical
                        counter[17] = counter[17] + 1;
                    }
                    if (protein1[i] == stop)
                    {
#pragma omp critical
                        counter[18] = counter[18] + 1;
                    }
                    if (protein1[i] == C)
                    {
#pragma omp critical
                        counter[19] = counter[19] + 1;
                    }
                    if (protein1[i] == W)
                    {
#pragma omp critical
                        counter[20] = counter[20] + 1;
                    }
                }
            }
        }
    }
    end = omp_get_wtime();
    for (int i = 0; i < 21; i++)
    {
        percent = (float)(((float)counter[i] / (float)lengthP) * 100);
        sum = sum + percent;
        cout << "Frequency of:  " << letters[i] << " is " << counter[i] << " Percentage is  " << percent << "%" << "\n";
    }

    cout << "Time" << end - start << "\n";
    cout << "Length of protein: " << lengthP << "\n";
    cout << "Total Frequency: " << sum << "%" << "\n";


}

int main()
{
    int option = 0;
    cout << "Please Choose One Of the following Options:" << "\n"; cout << "1) Sequential Mode" << "\n";
    cout << "2) Parallel For Mode" << "\n" << "3) Sections Mode" << "\n";
    cout << "Your Choice:  " << "\n"; cin >> option;
    if (option == 1)
    {
        cout << "The DNA Sequence Entered: " << "\n";
        read_seq();
        cout << "The Translated Protein: " << "\n";
        serialcode();
    }
    if (option == 2)
    {
        cout << "The DNA Sequence Entered: " << "\n";
        read_seq();
        cout << "The Translated Protein: " << "\n";
        paralelcode();
    }
    if (option == 3)
    {
        cout << "The DNA Sequence Entered: " << "\n";
        read_seq();
        cout << "The Translated Protein: " << "\n";
        prallelSectionCode();
    }

}


