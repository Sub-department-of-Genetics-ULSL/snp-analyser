const rawBaseUrl = import.meta.env.VITE_API_BASE_URL || '/api';
const baseUrl = rawBaseUrl.replace(/\/$/, '');

function toApiUrl(path) {
  return `${baseUrl}${path}`;
}

function toReportUrl(path) {
  if (!path) {
    return '';
  }

  if (/^https?:\/\//i.test(path)) {
    return path;
  }

  if (/^https?:\/\//i.test(baseUrl)) {
    return `${new URL(baseUrl).origin}${path}`;
  }

  return path;
}

async function parseJson(response, errorLabel) {
  if (!response.ok) {
    throw new Error(errorLabel);
  }
  return response.json();
}

/** The `detail` message FastAPI sends with an error response, if there is one. */
async function errorDetail(response, fallback) {
  try {
    const body = await response.json();
    return typeof body.detail === 'string' ? body.detail : fallback;
  } catch {
    return fallback;
  }
}

export const api = {
  async getOrganisms() {
    const response = await fetch(toApiUrl('/organisms'));
    return parseJson(response, 'Failed to load organisms');
  },

  async getGenes(latinName) {
    const response = await fetch(toApiUrl(`/organisms/${latinName}`));
    return parseJson(response, 'Failed to load genes');
  },

  async getGeneSequence(latinName, gene) {
    const response = await fetch(toApiUrl(`/organisms/${latinName}/${gene}`));
    return parseJson(response, 'Failed to load sequence');
  },

  async generateReport(organism, gene, mutations, predictMutatedStructure = false) {
    const response = await fetch(toApiUrl('/report'), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        organism,
        gene,
        mutations,
        predict_mutated_structure: predictMutatedStructure,
      })
    });

    if (!response.ok) {
      throw new Error('Failed to generate report');
    }

    return response.text();
  },

  toApiUrl,
  toReportUrl,

  async createReportJob(organism, gene, mutations, predictMutatedStructure = false) {
    const response = await fetch(toApiUrl('/report/jobs'), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        organism,
        gene,
        mutations,
        predict_mutated_structure: predictMutatedStructure,
      })
    });

    if (!response.ok) {
      throw new Error(await errorDetail(response, 'Failed to create report job'));
    }
    return response.json();
  },

  async describeVariants(organism, gene, mutations) {
    const response = await fetch(toApiUrl('/hgvs/describe'), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ organism, gene, mutations })
    });

    if (!response.ok) {
      throw new Error(await errorDetail(response, 'Failed to check the HGVS descriptions'));
    }
    return response.json();
  },

  async getReportJobs() {
    const response = await fetch(toApiUrl('/report/jobs'));
    return parseJson(response, 'Failed to load report queue');
  },

  async getReportJob(jobId) {
    const response = await fetch(toApiUrl(`/report/jobs/${jobId}`));
    return parseJson(response, 'Failed to load report job');
  }
};
