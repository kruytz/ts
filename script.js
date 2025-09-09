import VolleyballSimModule from './volleyball_sim.js';


const runButton = document.getElementById('runButton');
const numSimsInput = document.getElementById('numSims');
const statusDiv = document.getElementById('status');
const plotDiv = document.getElementById('plot');

runButton.addEventListener('click', runSimulation);

async function runSimulation() {
    runButton.disabled = true;
    runButton.textContent = 'Simulating...';
    statusDiv.textContent = 'Initializing WebAssembly module... Please wait.';
    plotDiv.innerHTML = '';

    try {
        const Module = await VolleyballSimModule();
        statusDiv.textContent = 'Module loaded. Starting simulations... This may take a moment.';

        const numSims = parseInt(numSimsInput.value, 10);
        console.log(`Running ${numSims} simulations...`);

        const startTime = performance.now();
        const csvData = Module.runSimulationsAndGetData(numSims);
        const endTime = performance.now();
        
        const duration = ((endTime - startTime) / 1000).toFixed(2);
        console.log(`Simulation finished in ${duration} seconds.`);
        
        const lines = csvData.trim().split('\n');
        lines.shift();
        const finalX = [];
        const finalZ = [];

        for (const line of lines) {
            const values = line.split(',');
            finalX.push(parseFloat(values[0]));
            finalZ.push(parseFloat(values[1]));
        }

        const validServesCount = finalX.length;
        statusDiv.innerHTML = `Found ${validServesCount} valid serves in ${duration} seconds. Calculating heatmap...`;

        if (validServesCount > 5) {
            createKDEHeatmap(finalX, finalZ, duration);
        } else {
            statusDiv.textContent = 'Not enough valid serves found to generate a heatmap. Try increasing the simulation count.';
        }

    } catch (error) {
        console.error("Error during simulation:", error);
        statusDiv.textContent = 'An error occurred. Check the console for details.';
    } finally {
        runButton.disabled = false;
        runButton.textContent = 'Run Simulation';
    }
}

function createKDEHeatmap(xData, zData, duration) {
    setTimeout(() => {
        const kde = window.science.stats.kde().sample(
            xData.map((x, i) => [x, zData[i]])
        );

        const n = 120;
        const m = 60;
        const xGrid = [], yGrid = [], zGrid = [];

        for (let i = 0; i < n; i++) xGrid.push(9 + (i * 9) / (n - 1));
        for (let j = 0; j < m; j++) yGrid.push(-4.5 + (j * 9) / (m - 1));

        for (let j = 0; j < m; j++) {
            const row = [];
            for (let i = 0; i < n; i++) {
                row.push(kde([xGrid[i], yGrid[j]]));
            }
            zGrid.push(row);
        }
        
        const trace = {
            x: xGrid,
            y: yGrid,
            z: zGrid,
            type: 'contour',
            colorscale: 'Inferno',
            contours: { coloring: 'heatmap' },
            hoverinfo: 'none'
        };

        const layout = {
            title: `Heatmap of ${xData.length} Valid Serves (High-Detail KDE)`,
            xaxis: { title: 'Horizontal Distance from Serve Line (m)', range: [-1, 19], scaleanchor: "y", scaleratio: 1 },
            yaxis: { title: 'Sideways Distance from Center (m)', range: [-5.5, 5.5] },
            shapes: [
                {type: 'rect', x0: 9, y0: -4.5, x1: 18, y1: 4.5, line: {color: 'black', width: 2}},
                {type: 'rect', x0: 0, y0: -4.5, x1: 9, y1: 4.5, line: {color: 'rgba(0,0,0,0.5)', width: 1}, fillcolor: 'rgba(230, 230, 230, 0.3)'},
                {type: 'line', x0: 9, y0: -4.5, x1: 9, y1: 4.5, line: {color: 'black', width: 3}},
                {type: 'line', x0: 12, y0: -4.5, x1: 12, y1: 4.5, line: {color: 'black', width: 1, dash: 'dash'}},
                {type: 'line', x0: 6, y0: -4.5, x1: 6, y1: 4.5, line: {color: 'black', width: 1, dash: 'dash'}},
            ],
            autosize: true,
            margin: { l: 60, r: 20, t: 40, b: 50 },
        };

        Plotly.newPlot('plot', [trace], layout, {responsive: true});
        
        statusDiv.innerHTML = `<strong>Simulation Complete!</strong> Found ${xData.length} valid serves in ${duration} seconds.`;
    }, 20);
}