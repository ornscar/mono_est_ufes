# Modelos de regressão quantílica com fração de cura: uma aplicação aos dados de COVID-19 grave na população materna

Trabalho de Conclusão de Curso "Modelos de regressão quantílica com fração de cura: uma aplicação aos dados de COVID-19 grave na população materna" orientado pela Profa. Dra. Agatha Sacramento Rodrigues e coorientado pelo Prof. Dr. Patrick Borges. A banca examinadora foi composta pela Profa. Dra. Rossana Pulcineli Vieira Francisco (FMUSP) e pelo Dr. Vinicius Fernando Calsavara (Cedars-Sinai Medical Center - Los Angeles).

### Resumo

Neste trabalho, abordamos modelos paramétricos de regressão quantílica para dados de sobrevivência com possibilidade de cura em que as distribuições são convenientemente reparametrizadas em termos do q-ésimo quantil ligado a covariáveis por intermédio de uma função logarítmica. Desenvolvemos o modelo de regressão quantílica de mistura padrão com distribuição Gompertz generalizada e, por meio de um estudo de simulação de Monte Carlo, comparamos com o modelo que considera a distribuição Gompertz generalizada em uma versão defeituosa (por isso, modelo defeituoso). Mostramos que os parâmetros do modelo defeituoso não ficam sobrecarregados ao estimar, simultaneamente, a fração de cura e os parâmetros do tempo de vida dos indivíduos que estão sujeitos à falha, pois seus resultados de simulação foram melhores em comparação ao modelo de mistura padrão quando os dados foram gerados pelo modelo defeituoso. Além disso, o modelo defeituoso e de mistura foram aplicados aos dados públicos de gestantes e puérperas de 10 a 55 anos internadas com Síndrome Respiratória Aguda Grave por COVID-19 que estão disponíveis no portal openDataSUS, do Ministério da Saúde. Portanto, este estudo ainda possibilitou analisar os efeitos de variáveis de caracterização, sintoma e comorbidade na fração de cura e em diferentes quantis dos tempos de sobrevivência dessas mulheres que estavam vivendo a gestação ou puerpério durante a pandemia de COVID-19 no Brasil.

### Script

Ajuste dos modelos de mistura padrão e defeituosos - `sadasdaaplicao/`

Apresentação - `apresentacao/`

### Software

R, versão 4.2.2, sob a IDE RStudio

### Sistema Operacional

Windows 10 Home Single Language, com Processador Intel Core i7 Octa-Core 1,8GHz e 16GB RAM
