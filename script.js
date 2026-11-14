const exampleMatrix = [[2,1,3,2,900],[2,4,1,3,1020],[-1,-1,1,0,20],[-1,0,0,2,0]];
const exampleNames = ['A','B','C','D'];
const tolerance = 1e-9;
let dimension = 4;
let selectedMethod = 'gauss';
let lastResult = null;
let lastMatrix = null;

const matrixEditor = document.getElementById('matrixEditor');
const sizeSelect = document.getElementById('sizeSelect');
const variableLabels = document.getElementById('variableLabels');
const errorBox = document.getElementById('errorBox');
const methodNames = {gauss:'Eliminación de Gauss',jordan:'Gauss-Jordan',cramer:'Regla de Cramer',inverse:'Matriz inversa'};

function variableName(index){return index<exampleNames.length?exampleNames[index]:`x${index+1}`}
function renderEditor(values){
  dimension = Number(sizeSelect.value);
  matrixEditor.innerHTML = '';
  matrixEditor.style.gridTemplateColumns = `repeat(${dimension}, minmax(42px, 1fr)) 12px minmax(58px, 1.15fr)`;
  for(let i=0;i<dimension;i++){
    const row=document.createElement('div');
    row.className='matrix-row';
    row.style.gridColumn=`1 / -1`;
    row.style.display='grid';
    row.style.gridTemplateColumns=`repeat(${dimension}, minmax(42px, 1fr)) 12px minmax(58px, 1.15fr)`;
    for(let j=0;j<=dimension;j++){
      if(j===dimension){const separator=document.createElement('span');separator.textContent='|';separator.style.textAlign='center';separator.style.color='#d4a477';row.appendChild(separator)}
      const input=document.createElement('input');
      input.type='number';input.step='any';input.className=`matrix-cell ${j===dimension?'rhs':''}`;input.dataset.row=i;input.dataset.col=j;
      input.setAttribute('aria-label',j===dimension?`Término independiente fila ${i+1}`:`Coeficiente ${variableName(j)} fila ${i+1}`);
      input.value=values&&values[i]&&values[i][j]!==undefined?values[i][j]:0;
      input.addEventListener('input',()=>{clearError();document.getElementById('inputStatus').textContent=`Sistema ${dimension} × ${dimension} · Datos editados · Pendiente de resolver`});
      row.appendChild(input);
    }
    matrixEditor.appendChild(row);
  }
  variableLabels.innerHTML=Array.from({length:dimension},(_,i)=>`<span>${variableName(i)}</span>`).join('');
  document.getElementById('inputStatus').textContent=`Sistema ${dimension} × ${dimension} · ${dimension} ecuaciones · ${dimension} incógnitas`;
}
function getMatrix(){
  const matrix=Array.from({length:dimension},()=>Array(dimension+1).fill(0));
  document.querySelectorAll('.matrix-cell').forEach(input=>{const value=input.value.trim();if(value==='')throw new Error('Todos los campos deben tener un valor numérico.');const number=Number(value);if(!Number.isFinite(number))throw new Error('Se encontró un valor que no es un número válido.');matrix[Number(input.dataset.row)][Number(input.dataset.col)]=number});
  if(matrix.every(row=>row.slice(0,dimension).every(v=>Math.abs(v)<tolerance))){throw new Error('La matriz de coeficientes no puede contener únicamente ceros.')}
  return matrix;
}
function clearError(){errorBox.hidden=true;errorBox.textContent=''}
function showError(message){errorBox.textContent=message;errorBox.hidden=false;errorBox.scrollIntoView({behavior:'smooth',block:'center'})}
function clone(matrix){return matrix.map(row=>row.slice())}
function fmt(value){if(Math.abs(value)<1e-10)return '0';if(Math.abs(value-Math.round(value))<1e-9)return String(Math.round(value));return Number(value.toFixed(6)).toString()}
function matrixHTML(matrix,augmented=true){const cols=augmented?matrix[0].length:matrix[0].length;return `<div class="matrix-display">${matrix.map(row=>`<div class="matrix-display-row">${row.map((v,j)=>`<span>${fmt(v)}</span>`).join('')}</div>`).join('')}</div>`}
function addStep(steps,title,description,matrix,formula){steps.push({title,description,matrix:matrix?clone(matrix):null,formula:formula||''})}
function rankOf(matrix,columns){const a=matrix.map(row=>row.slice(0,columns));let r=0;for(let c=0;c<columns&&r<a.length;c++){let p=r;for(let i=r+1;i<a.length;i++)if(Math.abs(a[i][c])>Math.abs(a[p][c]))p=i;if(Math.abs(a[p][c])<tolerance)continue;[a[p],a[r]]=[a[r],a[p]];const pivot=a[r][c];for(let j=c;j<columns;j++)a[r][j]/=pivot;for(let i=r+1;i<a.length;i++){const factor=a[i][c];for(let j=c;j<columns;j++)a[i][j]-=factor*a[r][j]}r++}return r}
function classify(matrix){const rankA=rankOf(matrix,dimension);const rankAug=rankOf(matrix,dimension+1);if(rankA<rankAug)return {type:'inconsistent',rankA,rankAug};if(rankA<dimension)return {type:'infinite',rankA,rankAug};return {type:'unique',rankA,rankAug}}
function gauss(matrix,steps){
  const a=clone(matrix),n=dimension;let pivotRow=0;
  addStep(steps,'Matriz aumentada inicial','Se organiza el sistema en la matriz aumentada [A | B].',a);
  for(let col=0;col<n&&pivotRow<n;col++){
    let pivot=pivotRow;for(let r=pivotRow+1;r<n;r++)if(Math.abs(a[r][col])>Math.abs(a[pivot][col]))pivot=r;
    if(Math.abs(a[pivot][col])<tolerance)continue;
    if(pivot!==pivotRow){[a[pivot],a[pivotRow]]=[a[pivotRow],a[pivot]];addStep(steps,`Intercambio F${pivotRow+1} ↔ F${pivot+1}`,`Se intercambian las filas para elegir un pivote de mayor valor absoluto y mejorar la estabilidad numérica.`,a)}
    const pivotValue=a[pivotRow][col];
    for(let r=pivotRow+1;r<n;r++){
      const factor=a[r][col]/pivotValue;if(Math.abs(factor)<tolerance)continue;
      for(let c=col;c<=n;c++)a[r][c]-=factor*a[pivotRow][c];
      a[r][col]=0;addStep(steps,`Operación F${r+1} ← F${r+1} − (${fmt(factor)})F${pivotRow+1}`,`Se elimina el elemento de la columna ${col+1} debajo del pivote.`,a);
    }
    pivotRow++;
  }
  const type=classify(a);
  if(type.type!=='unique')return {type,...type,steps,matrix:a};
  const x=Array(n).fill(0);
  for(let i=n-1;i>=0;i--){let sum=a[i][n];for(let j=i+1;j<n;j++)sum-=a[i][j]*x[j];x[i]=sum/a[i][i];addStep(steps,`Sustitución regresiva: ${variableName(i)}`,`Se despeja ${variableName(i)} usando las incógnitas calculadas en las filas anteriores.`,null,`${variableName(i)} = ${fmt(x[i])}`)}
  return {type:'unique',x,steps,matrix:a};
}
function jordan(matrix,steps){
  const a=clone(matrix),n=dimension;let pivotRow=0;
  addStep(steps,'Matriz aumentada inicial','Se parte de la matriz [A | B] y se transforma hasta obtener la identidad en el bloque de coeficientes.',a);
  for(let col=0;col<n&&pivotRow<n;col++){
    let pivot=pivotRow;for(let r=pivotRow+1;r<n;r++)if(Math.abs(a[r][col])>Math.abs(a[pivot][col]))pivot=r;
    if(Math.abs(a[pivot][col])<tolerance)continue;
    if(pivot!==pivotRow){[a[pivot],a[pivotRow]]=[a[pivotRow],a[pivot]];addStep(steps,`Intercambio F${pivotRow+1} ↔ F${pivot+1}`,'Se ubica un pivote no nulo en la fila de trabajo.',a)}
    const pv=a[pivotRow][col];for(let c=0;c<=n;c++)a[pivotRow][c]/=pv;
    addStep(steps,`Normalización F${pivotRow+1} ← F${pivotRow+1} / ${fmt(pv)}`,'Se convierte el pivote en uno dividiendo toda la fila por su valor.',a);
    for(let r=0;r<n;r++)if(r!==pivotRow){const factor=a[r][col];if(Math.abs(factor)<tolerance)continue;for(let c=0;c<=n;c++)a[r][c]-=factor*a[pivotRow][c];a[r][col]=0;addStep(steps,`Eliminación F${r+1} ← F${r+1} − (${fmt(factor)})F${pivotRow+1}`,'Se hacen ceros tanto por encima como por debajo del pivote.',a)}
    pivotRow++;
  }
  const type=classify(a);if(type.type!=='unique')return {type:'unique'===type.type?'unique':type.type,...type,steps,matrix:a};
  const x=a.map(row=>row[n]);addStep(steps,'Matriz escalonada reducida','El bloque izquierdo queda como la matriz identidad; la última columna contiene la solución.',a);
  return {type:'unique',x,steps,matrix:a};
}
function determinant(matrix){const n=matrix.length;if(n===1)return matrix[0][0];const a=clone(matrix);let det=1,sign=1;for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(a[r][c])>Math.abs(a[p][c]))p=r;if(Math.abs(a[p][c])<tolerance)return 0;if(p!==c){[a[p],a[c]]=[a[c],a[p]];sign*=-1}const pivot=a[c][c];det*=pivot;for(let r=c+1;r<n;r++){const factor=a[r][c]/pivot;for(let j=c+1;j<n;j++)a[r][j]-=factor*a[c][j]}}return det*sign}
function coefficientMatrix(matrix){return matrix.map(row=>row.slice(0,dimension))}
function cramer(matrix,steps){const A=coefficientMatrix(matrix),b=matrix.map(row=>row[dimension]),detA=determinant(A);addStep(steps,'Matriz de coeficientes y vector independiente','La regla de Cramer requiere una matriz cuadrada con determinante diferente de cero.',A);addStep(steps,'Determinante principal','Se calcula el determinante de la matriz de coeficientes.',null,`det(A) = ${fmt(detA)}`);if(Math.abs(detA)<tolerance){const c=classify(matrix);return {type:c.type,steps,matrix,detA}}const x=[];for(let j=0;j<dimension;j++){const Aj=clone(A);for(let i=0;i<dimension;i++)Aj[i][j]=b[i];const detAj=determinant(Aj);x[j]=detAj/detA;addStep(steps,`Determinante para ${variableName(j)}` ,`Se reemplaza la columna ${j+1} por el vector de términos independientes.`,Aj,`${variableName(j)} = det(A${j+1}) / det(A) = ${fmt(detAj)} / ${fmt(detA)} = ${fmt(x[j])}`)}return {type:'unique',x,steps,detA}}
function inverse(matrix,steps){const A=coefficientMatrix(matrix),n=dimension,b=matrix.map(row=>row[n]),detA=determinant(A);addStep(steps,'Matriz de coeficientes A','Se separa la matriz de coeficientes del vector de términos independientes.',A);addStep(steps,'Determinante de A','La matriz inversa existe únicamente si el determinante es distinto de cero.',null,`det(A) = ${fmt(detA)}`);if(Math.abs(detA)<tolerance){const c=classify(matrix);return {type:c.type,steps,matrix,detA}}const aug=A.map((row,i)=>[...row,...Array.from({length:n},(_,j)=>i===j?1:0)]);addStep(steps,'Matriz aumentada [A | I]','Se adjunta la matriz identidad a la matriz de coeficientes.',aug);for(let col=0;col<n;col++){let p=col;for(let r=col+1;r<n;r++)if(Math.abs(aug[r][col])>Math.abs(aug[p][col]))p=r;if(Math.abs(aug[p][col])<tolerance)return {type:'inconsistent',steps};if(p!==col){[aug[p],aug[col]]=[aug[col],aug[p]];addStep(steps,`Intercambio F${col+1} ↔ F${p+1}`,'Se intercambian filas para obtener un pivote no nulo.',aug)}const pv=aug[col][col];for(let j=0;j<2*n;j++)aug[col][j]/=pv;addStep(steps,`Normalización de F${col+1}`,'Se convierte el pivote en uno.',aug);for(let r=0;r<n;r++)if(r!==col){const f=aug[r][col];for(let j=0;j<2*n;j++)aug[r][j]-=f*aug[col][j];aug[r][col]=0;addStep(steps,`Eliminación de F${r+1}`,'Se obtiene un cero en la columna del pivote.',aug)}}const inv=aug.map(row=>row.slice(n));addStep(steps,'Matriz inversa A⁻¹','El bloque derecho de la matriz reducida corresponde a la inversa.',inv);const x=inv.map(row=>row.reduce((sum,v,j)=>sum+v*b[j],0));addStep(steps,'Producto A⁻¹ · B','Se multiplica la matriz inversa por el vector de términos independientes.',null,`X = A⁻¹ · B = [${x.map(fmt).join(', ')}]ᵀ`);return {type:'unique',x,steps,detA,inverse:inv}}
function residuals(matrix,x){return matrix.map((row,i)=>{const calculated=row.slice(0,dimension).reduce((sum,v,j)=>sum+v*x[j],0);return {equation:i+1,calculated,expected:row[dimension],error:Math.abs(calculated-row[dimension])}})}
function solve(method,matrix){const steps=[];let result;if(method==='gauss')result=gauss(matrix,steps);if(method==='jordan')result=jordan(matrix,steps);if(method==='cramer')result=cramer(matrix,steps);if(method==='inverse')result=inverse(matrix,steps);if(result.type==='unique'){result.residuals=residuals(matrix,result.x);result.verified=result.residuals.every(r=>r.error<1e-6)}return result}
function renderSteps(steps){return steps.map((step,index)=>`<article class="step-card"><div class="step-title"><span class="step-index">${String(index+1).padStart(2,'0')}</span><b>${step.title}</b></div><div class="step-body"><p>${step.description}</p>${step.matrix?matrixHTML(step.matrix):''}${step.formula?`<div class="inline-math">${step.formula}</div>`:''}</div></article>`).join('')}
function renderSolution(result,matrix){if(result.type!=='unique')return `<article class="step-card"><div class="step-body"><h3>${result.type==='infinite'?'Sistema compatible indeterminado':'Sistema incompatible'}</h3><p>${result.type==='infinite'?`El sistema tiene infinitas soluciones porque el rango de A (${result.rankA??rankOf(matrix,dimension)}) es menor que el número de incógnitas (${dimension}).`:'El sistema no tiene solución porque existe una contradicción entre las ecuaciones.'}</p>${result.type==='infinite'?`<p>Para describir el conjunto solución se requieren variables libres y parámetros. Los métodos Cramer e inversa no producen una solución única en este caso.</p>`:''}</div></article>`;
 const vars=result.x.map((v,i)=>`<div class="solution-item"><span>${variableName(i)}</span><b>${fmt(v)}</b></div>`).join('');const rows=result.residuals.map(r=>`<tr><td>Ecuación ${r.equation}</td><td>${fmt(r.calculated)}</td><td>${fmt(r.expected)}</td><td class="${r.error<1e-6?'valid':'invalid'}">${r.error<1e-6?'Correcta':'Revisar'} · ${fmt(r.error)}</td></tr>`).join('');return `<h3>Valores de las incógnitas</h3><div class="solution-grid">${vars}</div><div class="step-card"><div class="step-title"><span class="step-index">✓</span><b>Comprobación por sustitución</b></div><div class="step-body"><p>Se sustituyen los valores obtenidos en cada ecuación original y se compara el lado izquierdo con el término independiente.</p><table class="validation-table"><thead><tr><th>Ecuación</th><th>Resultado calculado</th><th>Valor esperado</th><th>Verificación</th></tr></thead><tbody>${rows}</tbody></table><p style="margin-top:13px"><b class="${result.verified?'valid':'invalid'}">${result.verified?'La solución satisface todas las ecuaciones.':'La solución requiere revisión numérica.'}</b></p></div></div>${matrix.length===4&&dimension===4?`<div class="step-card"><div class="step-title"><span class="step-index">↗</span><b>Interpretación en el contexto de producción</b></div><div class="step-body"><p>De acuerdo con el modelo, la cantidad calculada de cada compuesto es:</p><p>${result.x.map((v,i)=>`${variableName(i)}: ${fmt(v)} unidades`).join(' · ')}.</p><p>Las cantidades satisfacen las disponibilidades de MP1 y MP2 y las dos condiciones adicionales establecidas en el planteamiento.</p></div></div>`:''}`}
function renderComparison(matrix){const methods=['gauss','jordan','cramer','inverse'];return `<div class="compare-grid">${methods.map(method=>{const r=solve(method,matrix);return `<article class="compare-item"><h4>${methodNames[method]}</h4><p>${r.type==='unique'?`Solución: ${r.x.map((v,i)=>`${variableName(i)} = ${fmt(v)}`).join(', ')}. Verificación: ${r.verified?'correcta':'pendiente'}.`:`${r.type==='infinite'?'Sistema con infinitas soluciones.':'Sistema incompatible o sin solución única.'} ${['cramer','inverse'].includes(method)?'No se puede obtener una solución única mediante este método.':''}`}</p><p>Operaciones registradas: ${r.steps.length}</p></article>`}).join('')}</div><p style="font-size:10px;color:#839087;margin-top:14px">Los métodos se ejecutan sobre la misma matriz. Gauss y Gauss-Jordan permiten analizar sistemas no únicos; Cramer e inversa requieren una matriz cuadrada no singular para obtener una solución única.</p>`}
function displayResult(method,matrix,result){lastResult=result;lastMatrix=clone(matrix);document.getElementById('emptyState').hidden=true;document.getElementById('resultsContent').hidden=false;const title=result.type==='unique'?'Solución única encontrada':result.type==='infinite'?'Sistema con infinitas soluciones':'Sistema incompatible';const detail=result.type==='unique'?`El método ${methodNames[method]} encontró una solución y se verificó mediante sustitución.`:result.type==='infinite'?'El sistema es compatible indeterminado; no existe una solución única.':'El sistema no tiene solución que satisfaga simultáneamente todas las ecuaciones.';document.getElementById('resultSummary').innerHTML=`<div><h3>${title}</h3><p>${detail}</p></div><span class="result-badge">${methodNames[method]}</span>`;document.getElementById('stepsPanel').innerHTML=renderSteps(result.steps);document.getElementById('solutionPanel').innerHTML=renderSolution(result,matrix);document.getElementById('comparisonPanel').innerHTML=renderComparison(matrix);setTab('steps');document.getElementById('resultados').scrollIntoView({behavior:'smooth',block:'start'})}
function setTab(tab){document.querySelectorAll('.result-tab').forEach(button=>button.classList.toggle('active',button.dataset.tab===tab));['steps','solution','comparison'].forEach(name=>document.getElementById(`${name}Panel`).hidden=name!==tab)}
document.querySelectorAll('.method-option').forEach(button=>button.addEventListener('click',()=>{selectedMethod=button.dataset.method;document.querySelectorAll('.method-option').forEach(item=>item.classList.toggle('selected',item===button));clearError()}));
document.querySelectorAll('.result-tab').forEach(button=>button.addEventListener('click',()=>setTab(button.dataset.tab)));
sizeSelect.addEventListener('change',()=>{const next=Number(sizeSelect.value);const old=getMatrixSafe();const resized=Array.from({length:next},(_,i)=>Array.from({length:next+1},(_,j)=>old[i]?.[j]??0));renderEditor(resized);clearError()});
function getMatrixSafe(){try{return getMatrix()}catch{return Array.from({length:dimension},()=>Array(dimension+1).fill(0))}}
document.getElementById('loadExample').addEventListener('click',()=>{sizeSelect.value='4';renderEditor(exampleMatrix);clearError();document.getElementById('inputStatus').textContent='Ejemplo cargado · Sistema 4 × 4'});
document.getElementById('resetButton').addEventListener('click',()=>{renderEditor(Array.from({length:dimension},()=>Array(dimension+1).fill(0)));clearError();document.getElementById('resultsContent').hidden=true;document.getElementById('emptyState').hidden=false});
document.getElementById('solveButton').addEventListener('click',()=>{try{clearError();const matrix=getMatrix();const result=solve(selectedMethod,matrix);displayResult(selectedMethod,matrix,result)}catch(error){showError(error.message)}});
renderEditor(exampleMatrix);
