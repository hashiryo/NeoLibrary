/**
 * 自作 SSG ビルドスクリプト
 * VitePress を使わず、1ページずつ処理して OOM を回避する
 *
 * 実行: npx tsx scripts/build-site.ts
 */
import fs from "fs";
import path from "path";
import MarkdownIt from "markdown-it";
import { katex } from "@mdit/plugin-katex";
import matter from "gray-matter";
import { createHighlighter, type Highlighter } from "shiki";
import katexLib from "katex";
import {
  buildDependencyGraph,
  type DependencyGraph,
} from "./lib/dependency-graph";
import { bundleCpp } from "./lib/bundle";

// ============================================================
// 定数
// ============================================================

const ROOT = path.resolve(__dirname, "..");
const SRC_DIR = path.join(ROOT, "neo");
const MD_DIR = path.join(ROOT, "md");
const SITE_DIR = path.join(ROOT, "site");
const OUT_DIR = path.join(SITE_DIR, "NeoLibrary");
const BASE_PATH = "/NeoLibrary";
// procon-judge のサイト。ヘッダごとの逆引き JSON (data/headers/<path>.json) を
// 表示時に読む。ビルド時に焼き込まないのは順番の都合で、neo を push すると
// こちらのサイトが先に建ち、judge が測るのはそのあとになるため。焼き込むと
// 定常状態で常に「参考」と出る。
const JUDGE_SITE = "https://hashiryo.github.io/procon-judge";

// ============================================================
// markdown-it 初期化
// ============================================================

let highlighter: Highlighter;

async function initMarkdown() {
  highlighter = await createHighlighter({
    themes: ["github-light", "github-dark"],
    langs: ["cpp"],
  });

  const md = MarkdownIt({ html: true });
  md.use(katex);

  // shiki によるコードブロックハイライト
  const defaultFence = md.renderer.rules.fence!;
  md.renderer.rules.fence = (tokens, idx, options, env, self) => {
    const token = tokens[idx];
    const lang = token.info.trim();
    if (lang === "cpp" || lang === "c++") {
      return highlighter.codeToHtml(token.content, {
        lang: "cpp",
        themes: { light: "github-light", dark: "github-dark" },
      });
    }
    return defaultFence(tokens, idx, options, env, self);
  };

  return md;
}

// ============================================================
// HTML テンプレート
// ============================================================

/** テキスト内の $...$ をKaTeX HTMLに変換 */
function renderInlineKatex(text: string): string {
  return text.replace(/\$([^$]+)\$/g, (_, expr) => {
    try {
      return katexLib.renderToString(expr, { throwOnError: false });
    } catch {
      return escapeHtml(expr);
    }
  });
}

function escapeHtml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/**
 * レンダリング済みHTMLのリンク・画像パスを書き換える。
 * - 絶対パス /img/... → BASE_PATH/img/...
 * - 相対 .md リンク → 対応する .html へ
 * - 相対 .hpp リンク → サイト上の対応ページへ
 * mdDir: md ファイルのあるディレクトリ (例: "string")
 */
function rewriteLinks(html: string, mdDir: string): string {
  // href="..." と src="..." の両方を書き換える
  return html.replace(
    /((?:href|src)=")([^"]+)(")/g,
    (_match, pre, url, post) => {
      // 外部URLはスキップ
      if (/^https?:\/\//.test(url)) return pre + url + post;
      // # のみのアンカーはスキップ
      if (url.startsWith("#")) return pre + url + post;

      // 絶対パス → BASE_PATH を付与
      if (url.startsWith("/")) {
        return pre + BASE_PATH + url + post;
      }

      // 相対パスを md ディレクトリ基準で解決
      const resolved = path.posix.normalize(path.posix.join(mdDir, url));

      // .hpp → サイト上の対応ページ
      if (resolved.endsWith(".hpp")) {
        // neo/ prefix を除去して .html に変換
        const pagePath = resolved
          .replace(/^(?:\.\.\/)*neo\//, "")
          .replace(/\.hpp$/, ".html");
        return pre + `${BASE_PATH}/${pagePath}` + post;
      }

      // .md → .html に変換
      if (resolved.endsWith(".md")) {
        const pagePath = resolved.replace(/\.md$/, ".html");
        return pre + `${BASE_PATH}/${pagePath}` + post;
      }

      return pre + url + post;
    },
  );
}

/**
 * ヘッダの状態の印。中身は表示時に procon-judge の data/headers/index.json を読んで
 * 色を付ける (renderPage のスクリプト)。ビルド時に焼き込まないのは、NeoLibrary の
 * サイトが先に建って judge がそのあと測る順番のため。焼くと定常状態で常に古く出る。
 */
function judgeDot(hppPath: string): string {
  return `<span data-pagefind-ignore class="dot dot-gray" data-judge-dot="${escapeHtml(hppPath)}" title="procon-judge の記録を読んでいます">●</span>`;
}

// ============================================================
// レイアウト
// ============================================================

function renderPage(title: string, content: string, sidebar: string): string {
  return `<!DOCTYPE html>
<html lang="ja">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>${escapeHtml(title)} | Hashiryo's NeoLibrary</title>
  <link rel="stylesheet" href="${BASE_PATH}/assets/style.css">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
  <link rel="stylesheet" href="${BASE_PATH}/pagefind/pagefind-ui.css">
</head>
<body>
  <nav class="top-nav">
    <a href="${BASE_PATH}/" class="site-title">Hashiryo's NeoLibrary</a>
    <div id="search"></div>
    <a href="https://github.com/hashiryo/NeoLibrary">GitHub</a>
  </nav>
  <div class="layout">
    <aside class="sidebar">
      ${sidebar}
    </aside>
    <main class="content" data-pagefind-body>
      ${content}
    </main>
  </div>
  <script src="${BASE_PATH}/pagefind/pagefind-ui.js"></script>
  <script>
    new PagefindUI({ element: "#search", showSubResults: true, showImages: false });
    document.querySelectorAll('pre.shiki').forEach(pre => {
      const wrapper = document.createElement('div');
      wrapper.className = 'shiki-wrapper';
      pre.parentNode.insertBefore(wrapper, pre);
      wrapper.appendChild(pre);
      const btn = document.createElement('button');
      btn.className = 'copy-btn';
      btn.textContent = 'Copy';
      btn.addEventListener('click', () => {
        navigator.clipboard.writeText(pre.textContent).then(() => {
          btn.textContent = 'Copied!';
          setTimeout(() => btn.textContent = 'Copy', 2000);
        });
      });
      wrapper.appendChild(btn);
    });
    // Bundle toggle
    const toggle = document.querySelector('.code-toggle');
    if (toggle) {
      const original = document.getElementById('code-original');
      const bundled = document.getElementById('code-bundled');
      const bundleBtn = document.createElement('button');
      bundleBtn.className = 'bundle-btn';
      bundleBtn.textContent = 'Bundle';
      bundleBtn.addEventListener('click', () => {
        const showingOriginal = original.style.display !== 'none';
        original.style.display = showingOriginal ? 'none' : '';
        bundled.style.display = showingOriginal ? '' : 'none';
        bundleBtn.textContent = showingOriginal ? 'Original' : 'Bundle';
      });
      toggle.insertBefore(bundleBtn, toggle.firstChild);
    }
  </script>
  <script>
    // procon-judge の記録を表示時に読む。上のスクリプトとは分けておく。検索 UI が
    // 読めない回に巻き込まれないため。
    const JUDGE = '${JUDGE_SITE}';

    // ヘッダの状態の印。data/headers/index.json (ヘッダごとの要約) を 1 回読んで、
    // 印ごとに色と説明を付ける。読めなければ灰色のまま。
    const dots = document.querySelectorAll('[data-judge-dot]');
    if (dots.length) {
      fetch(JUDGE + '/data/headers/index.json')
        .then((r) => (r.ok ? r.json() : null))
        .then((data) => {
          if (!data || (data.schema ?? 1) !== 1 || !data.headers) return;
          dots.forEach((dot) => {
            const h = data.headers[dot.dataset.judgeDot];
            if (!h) {
              dot.title = 'このヘッダを使う提出は procon-judge に無い';
              return;
            }
            let ac = 0, failing = 0, stale = 0;
            for (const e of h.envs) { ac += e.ac; failing += e.failing; stale += e.stale; }
            const verified = h.envs.filter((e) => e.verified).length;
            let cls = 'dot-gray';
            if (h.verified) cls = 'dot-ac';
            else if (failing > 0 && ac === 0) cls = 'dot-fail';
            else if (ac > 0 || failing > 0) cls = 'dot-warn';
            // 実行時の分岐を持つヘッダは、Codeforces と同じ命令の CPU でも代わりの経路を
            // 確かめている (procon-judge の pj fallback)。通らない提出があれば緑にしない。
            const fb = h.fallback;
            if (fb && fb.failing > 0 && cls === 'dot-ac') cls = 'dot-warn';
            dot.className = 'dot ' + cls;
            let title;
            if (h.verified) title = '全環境で現行の AC';
            else if (ac === 0 && failing === 0) title = '今のコードでは未計測です';
            else title = verified + ' / ' + h.envs.length + ' 環境で現行の AC';
            if (failing) title += '、失敗 ' + failing;
            if (stale) title += '、参考 ' + stale;
            if (h.compile_only) title += ' (コンパイルのみ)';
            if (fb) {
              title += fb.failing
                ? '。Codeforces 相当の CPU で失敗 ' + fb.failing + ' (' + fb.checked + ' 本中)'
                : '。Codeforces 相当の CPU でも AC (' + fb.ac + ' 本)';
              if (fb.skipped) title += '、確かめられなかった ' + fb.skipped;
            }
            dot.title = title;
          });
        })
        .catch(() => {});
    }

    // Submissions (procon-judge)。ヘッダごとの逆引き JSON を表示時に読む。
    const judge = document.querySelector('.judge-section');
    if (judge) {
      const el = (tag, text, cls) => {
        const node = document.createElement(tag);
        if (text != null) node.textContent = text;
        if (cls) node.className = cls;
        return node;
      };
      const statusClass = (s) =>
        s === 'AC' ? 'status-ac'
        : s === 'WA' || s === 'RE' ? 'status-fail'
        : s === 'TLE' || s === 'MLE' ? 'status-warn'
        : 'status-gray';
      const table = (data, rows, site) => {
        const envs = data.environments || [];
        const t = el('table', null, 'verify-matrix judge-table');
        const head = el('tr');
        head.append(el('th', 'Problem'), el('th', 'Submission'));
        for (const env of envs) head.append(el('th', env));
        const thead = el('thead');
        thead.append(head);
        t.append(thead);
        const body = el('tbody');
        for (const s of rows) {
          const tr = el('tr');
          const problem = el('td');
          const pa = el('a', s.title || s.problem);
          pa.href = site + s.problem_page;
          problem.append(pa);
          const submission = el('td');
          const name = s.submission.startsWith('submissions/') ? s.submission.slice('submissions/'.length) : s.submission;
          const sa = el('a', name);
          sa.href = site + s.page;
          submission.append(sa);
          // 判定サイトのテストデータでない問題。AC の意味が違うので印を付ける。
          if (s.official === false) {
            const labels = { none: 'コンパイルのみ', local: '自作', manual: '手動取り込み' };
            // none は「走らせない」と「提出が自分で検証する」の 2 通りあるので、比べ方で分ける
            // (procon-judge のサイトの「無し (自己検証)」と同じ)。
            const selfCheck = s.testdata === 'none' && s.compare === 'exit_code';
            const chip = el('span', selfCheck ? '自己検証' : labels[s.testdata] || s.testdata, 'judge-chip');
            chip.title = selfCheck
              ? 'テストケースが無く、提出が自分で持っている入出力で検証します。終了コード 0 で走り切れば AC です'
              : 'テストケースは判定サイトのものではありません';
            submission.append(chip);
          }
          // Codeforces と同じ命令の CPU (QEMU) で、このヘッダの実行時の分岐の代わりの経路を
          // 走らせた結果 (procon-judge の pj fallback)。印を持つヘッダの、対象の提出にだけある。
          if (s.fallback) {
            const skipped = s.fallback === 'SKIP';
            const ok = s.fallback === 'AC' || skipped;
            const chip = el('span', 'CF 相当 ' + (skipped ? '未確認' : s.fallback), ok ? 'judge-chip' : 'judge-chip judge-chip-bad');
            chip.title = skipped
              ? 'テストデータを用意できず、Codeforces 相当の CPU では確かめられませんでした'
              : 'Codeforces と同じ命令の CPU (vpclmulqdq、GFNI、AVX-512 なし) で代わりの経路を走らせた結果: ' + s.fallback;
            submission.append(chip);
          }
          tr.append(problem, submission);
          const byEnv = Object.fromEntries((s.envs || []).map((e) => [e.env, e]));
          for (const env of envs) {
            const e = byEnv[env];
            if (!e || e.status == null) {
              tr.append(el('td', '-', 'status-gray'));
              continue;
            }
            const td = el('td', e.status, statusClass(e.status));
            if (e.current === false) {
              td.classList.add('judge-stale');
              td.title = '測ってからヘッダが変わっています。次の計測で入れ替わります';
              td.append(el('span', '参考', 'judge-chip'));
            }
            tr.append(td);
          }
          body.append(tr);
        }
        t.append(body);
        const wrap = el('div', null, 'table-wrapper');
        wrap.append(t);
        return wrap;
      };
      fetch(JUDGE + '/data/headers/' + encodeURI(judge.dataset.judgeHeader) + '.json')
        .then((r) => (r.ok ? r.json() : null))
        .then((data) => {
          if (!data || !data.submissions || data.submissions.length === 0) return;
          // JSON の契約の版。judge が形を変えて版を上げたら、知らない形を読まずに節を隠す。
          if ((data.schema ?? 1) !== 1) return;
          const site = data.site || JUDGE + '/';
          const direct = data.submissions.filter((s) => s.direct);
          const via = data.submissions.filter((s) => !s.direct);
          const out = judge.querySelector('.judge-body');
          if (direct.length) out.append(table(data, direct, site));
          if (via.length) {
            const details = el('details');
            details.append(el('summary', 'Indirect (' + via.length + ')'));
            details.append(table(data, via, site));
            out.append(details);
          }
          judge.hidden = false;
        })
        .catch(() => {});
    }
  </script>
</body>
</html>`;
}

// ============================================================
// サイドバー生成
// ============================================================

function readFrontmatter(mdPath: string): { title?: string; order?: number } {
  if (!fs.existsSync(mdPath)) return {};
  const content = fs.readFileSync(mdPath, "utf-8");
  const titleMatch = content.match(
    /^---\s*\n[\s\S]*?title:\s*(.+)\n[\s\S]*?---/,
  );
  const orderMatch = content.match(
    /^---\s*\n[\s\S]*?order:\s*(\d+)\n[\s\S]*?---/,
  );
  return {
    title: titleMatch ? titleMatch[1].trim() : undefined,
    order: orderMatch ? parseInt(orderMatch[1]) : undefined,
  };
}

function generateSidebar(): string {
  const categories = fs
    .readdirSync(SRC_DIR, { withFileTypes: true })
    .filter((e) => e.isDirectory());

  interface SidebarItem {
    text: string;
    icon?: string;
    link?: string;
    items?: SidebarItem[];
    order: number;
    dirName?: string;
  }

  function buildItems(
    srcPath: string,
    mdPath: string,
    linkPrefix: string,
    hppPrefix: string,
  ): SidebarItem[] {
    const entries = fs.readdirSync(srcPath, { withFileTypes: true });
    const items: SidebarItem[] = [];

    for (const entry of entries) {
      if (!entry.isFile() || !entry.name.endsWith(".hpp")) continue;
      const name = entry.name.replace(/\.hpp$/, "");
      const fm = readFrontmatter(path.join(mdPath, name + ".md"));
      const hppPath = `neo/${hppPrefix}${entry.name}`;
      const icon = judgeDot(hppPath);
      items.push({
        text: fm.title || name,
        icon,
        link: `${BASE_PATH}${linkPrefix}/${name}.html`,
        order: fm.order ?? 999,
      });
    }

    for (const entry of entries) {
      if (!entry.isDirectory()) continue;
      const subItems = buildItems(
        path.join(srcPath, entry.name),
        path.join(mdPath, entry.name),
        `${linkPrefix}/${entry.name}`,
        `${hppPrefix}${entry.name}/`,
      );
      if (subItems.length > 0) {
        const indexFm = readFrontmatter(
          path.join(mdPath, entry.name, "_index.md"),
        );
        items.push({
          text: indexFm.title || entry.name,
          dirName:
            indexFm.title && indexFm.title !== entry.name
              ? entry.name
              : undefined,
          items: subItems,
          order: indexFm.order ?? 999,
        });
      }
    }

    items.sort((a, b) =>
      a.order !== b.order ? a.order - b.order : a.text.localeCompare(b.text),
    );
    return items;
  }

  function renderItems(items: SidebarItem[]): string {
    let html = "<ul>";
    for (const item of items) {
      if (item.link) {
        const icon = item.icon ? `${item.icon} ` : "";
        html += `<li><a href="${item.link}">${icon}${renderInlineKatex(escapeHtml(item.text))}</a></li>`;
      } else if (item.items) {
        const dirLabel = item.dirName
          ? `<span class="sidebar-dir">${escapeHtml(item.dirName)}</span>`
          : "";
        html += `<li><details><summary>${renderInlineKatex(escapeHtml(item.text))}${dirLabel}</summary>${renderItems(item.items)}</details></li>`;
      }
    }
    html += "</ul>";
    return html;
  }

  const allItems: SidebarItem[] = [];
  for (const cat of categories) {
    const items = buildItems(
      path.join(SRC_DIR, cat.name),
      path.join(MD_DIR, cat.name),
      `/${cat.name}`,
      `${cat.name}/`,
    );
    if (items.length > 0) {
      const indexFm = readFrontmatter(path.join(MD_DIR, cat.name, "_index.md"));
      allItems.push({
        text: indexFm.title || cat.name,
        dirName:
          indexFm.title && indexFm.title !== cat.name ? cat.name : undefined,
        items,
        order: indexFm.order ?? 999,
      });
    }
  }
  allItems.sort((a, b) =>
    a.order !== b.order ? a.order - b.order : a.text.localeCompare(b.text),
  );

  return renderItems(allItems);
}

// ============================================================
// hpp ドキュメントページ生成
// ============================================================

function generateHppPage(
  md: MarkdownIt,
  hppPath: string,
  mdPath: string,
  sidebar: string,
  depGraph: DependencyGraph,
): void {
  const raw = fs.readFileSync(mdPath, "utf-8");
  const { data: fm, content: mdContent } = matter(raw);

  const icon = judgeDot(hppPath);
  const title =
    fm.title ||
    hppPath
      .split("/")
      .pop()
      ?.replace(/\.hpp$/, "") ||
    hppPath;

  // md コンテンツをレンダリング
  let body = `<h1>${icon} ${renderInlineKatex(escapeHtml(title))}</h1>\n`;

  // ソースコード (元 + バンドル切り替え)
  const sourcePath = path.join(ROOT, hppPath);
  if (fs.existsSync(sourcePath)) {
    const source = fs.readFileSync(sourcePath, "utf-8");
    const githubUrl = `https://github.com/hashiryo/NeoLibrary/blob/main/${hppPath}`;
    body += `<h2>Code</h2>\n`;
    body += `<details><summary>${escapeHtml(hppPath)}</summary>\n`;
    body += `<p><a href="${githubUrl}">View on GitHub</a></p>\n`;

    // 依存がある場合のみバンドルボタンを表示
    const hasDeps = (depGraph.dependsOn[hppPath] || []).length > 0;
    if (hasDeps) {
      const bundled = bundleCpp(path.resolve(ROOT, hppPath), {
        includeDirs: [ROOT],
        commentMarkers: true,
      });
      body += `<div class="code-toggle">`;
      body += `<div class="code-view active" id="code-original">`;
      body += highlighter.codeToHtml(source, {
        lang: "cpp",
        themes: { light: "github-light", dark: "github-dark" },
      });
      body += `</div>`;
      body += `<div class="code-view" id="code-bundled" style="display:none">`;
      body += highlighter.codeToHtml(bundled, {
        lang: "cpp",
        themes: { light: "github-light", dark: "github-dark" },
      });
      body += `</div>`;
      body += `</div>`;
    } else {
      body += highlighter.codeToHtml(source, {
        lang: "cpp",
        themes: { light: "github-light", dark: "github-dark" },
      });
    }
    body += `</details>\n`;
  }

  // ユーザーの md コンテンツ
  if (mdContent.trim()) {
    const mdRelDir = path.posix.dirname(path.relative(MD_DIR, mdPath));
    body += rewriteLinks(md.render(mdContent), mdRelDir);
  }

  // 正しさの証拠は procon-judge の記録だけ。Library で 2026-09-22 に verify を畳んだ形を
  // そのまま引き継ぐ (algo-notes の「Library の verify を畳む設計」)。

  // Submissions (procon-judge)
  // 中身は renderPage のスクリプトが表示時に judge から読んで埋める。JSON が
  // 無い (どの提出も使っていない) ヘッダでは節ごと隠したままにする。
  body += `<section class="judge-section" hidden data-judge-header="${escapeHtml(hppPath)}">\n`;
  body += "<h2>Submissions</h2>\n";
  body += `<p class="judge-note">procon-judge で計測した、このヘッダを使う提出。環境ごとの状態と、今のヘッダで測った記録かどうか。参考は測ってからヘッダが変わったもので、次の計測で入れ替わる。 <a href="${JUDGE_SITE}/">procon-judge</a></p>\n`;
  body += '<div class="judge-body"></div>\n</section>\n';

  function renderDepItem(hppRelPath: string): string {
    const icon = judgeDot(hppRelPath);
    const title =
      readFrontmatter(
        path.join(
          MD_DIR,
          hppRelPath.replace(/^neo\//, "").replace(/\.hpp$/, ".md"),
        ),
      ).title ||
      hppRelPath
        .split("/")
        .pop()
        ?.replace(/\.hpp$/, "") ||
      hppRelPath;
    const link = `${BASE_PATH}/${hppRelPath.replace(/^neo\//, "").replace(/\.hpp$/, ".html")}`;
    return `<li>${icon} <a href="${link}">${renderInlineKatex(escapeHtml(title))}</a><span class="dep-path">${escapeHtml(hppRelPath)}</span></li>\n`;
  }

  // Depends on
  const directDeps = new Set(depGraph.dependsOn[hppPath] || []);
  const allDeps = depGraph.transitiveDeps[hppPath] || new Set<string>();
  const indirectDeps = [...allDeps].filter((d) => !directDeps.has(d)).sort();
  if (directDeps.size > 0 || indirectDeps.length > 0) {
    body += "<h2>Depends on</h2>\n";
    if (directDeps.size > 0) {
      body += '<h3>Direct</h3>\n<ul class="dep-list">\n';
      for (const dep of [...directDeps].sort()) {
        body += renderDepItem(dep);
      }
      body += "</ul>\n";
    }
    if (indirectDeps.length > 0) {
      body += `<details><summary>Indirect (${indirectDeps.length})</summary>\n<ul class="dep-list">\n`;
      for (const dep of indirectDeps) {
        body += renderDepItem(dep);
      }
      body += "</ul>\n</details>\n";
    }
  }

  // Required by
  const directReqBy = new Set(depGraph.requiredBy[hppPath] || []);
  const allReqBy = depGraph.transitiveRequiredBy[hppPath] || new Set<string>();
  const indirectReqBy = [...allReqBy].filter((r) => !directReqBy.has(r)).sort();
  if (directReqBy.size > 0 || indirectReqBy.length > 0) {
    body += "<h2>Required by</h2>\n";
    if (directReqBy.size > 0) {
      body += '<h3>Direct</h3>\n<ul class="dep-list">\n';
      for (const req of [...directReqBy].sort()) {
        body += renderDepItem(req);
      }
      body += "</ul>\n";
    }
    if (indirectReqBy.length > 0) {
      body += `<details><summary>Indirect (${indirectReqBy.length})</summary>\n<ul class="dep-list">\n`;
      for (const req of indirectReqBy) {
        body += renderDepItem(req);
      }
      body += "</ul>\n</details>\n";
    }
  }

  // 出力
  const outRelPath = hppPath.replace(/^neo\//, "").replace(/\.hpp$/, ".html");
  const outPath = path.join(OUT_DIR, outRelPath);
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.writeFileSync(outPath, renderPage(title, body, sidebar));
}

// ============================================================
// ホームページ生成
// ============================================================

function generateHomePage(sidebar: string): void {
  const body = `
<div class="hero">
  <h1>Hashiryo's NeoLibrary</h1>
  <p>競技プログラミング用C++ライブラリ</p>
  <div class="hero-actions">
    <a href="https://github.com/hashiryo/NeoLibrary" class="btn">GitHub</a>
  </div>
</div>`;
  const outPath = path.join(OUT_DIR, "index.html");
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.writeFileSync(outPath, renderPage("Hashiryo's NeoLibrary", body, sidebar));
}

// ============================================================
// メイン
// ============================================================

async function main() {
  console.time("Total build");

  // 初期化
  const md = await initMarkdown();
  const depGraph = buildDependencyGraph();
  const sidebar = generateSidebar();

  // 出力ディレクトリを準備
  if (fs.existsSync(SITE_DIR)) fs.rmSync(SITE_DIR, { recursive: true });
  fs.mkdirSync(OUT_DIR, { recursive: true });

  // 静的アセットをコピー
  const publicDir = path.join(MD_DIR, "public");
  if (fs.existsSync(publicDir)) {
    function copyDir(src: string, dest: string) {
      fs.mkdirSync(dest, { recursive: true });
      for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
        const s = path.join(src, entry.name);
        const d = path.join(dest, entry.name);
        if (entry.isDirectory()) copyDir(s, d);
        else fs.copyFileSync(s, d);
      }
    }
    copyDir(publicDir, OUT_DIR);
  }

  // CSS を生成
  fs.mkdirSync(path.join(OUT_DIR, "assets"), { recursive: true });
  fs.writeFileSync(path.join(OUT_DIR, "assets", "style.css"), CSS);

  // ホームページ
  generateHomePage(sidebar);
  console.log("Generated: index.html");

  // hpp ドキュメントページ (1ページずつ処理)
  let hppCount = 0;
  function scanHpp(dir: string) {
    for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        scanHpp(full);
        continue;
      }
      if (!entry.name.endsWith(".hpp")) continue;

      const hppPath = path.relative(ROOT, full); // "neo/algebra/GF2p64.hpp"
      const mdRelPath = hppPath
        .replace(/^neo\//, "")
        .replace(/\.hpp$/, ".md");
      const mdPath = path.join(MD_DIR, mdRelPath);

      // md がなければスタブとして生成
      const mdExists = fs.existsSync(mdPath);
      if (!mdExists) {
        const name = entry.name.replace(/\.hpp$/, "");
        const stubMd = path.join(MD_DIR, mdRelPath);
        fs.mkdirSync(path.dirname(stubMd), { recursive: true });
        fs.writeFileSync(
          stubMd,
          `---\ntitle: ${name}\ndocumentation_of: ../../${hppPath}\n---\n`,
        );
      }

      generateHppPage(
        md,
        hppPath,
        mdPath || path.join(MD_DIR, mdRelPath),
        sidebar,
        depGraph,
      );
      hppCount++;
    }
  }
  scanHpp(SRC_DIR);
  console.log(`Generated: ${hppCount} hpp documentation pages`);

  console.timeEnd("Total build");
}

// ============================================================
// CSS
// ============================================================

const CSS = `
:root {
  --c-bg: #fff;
  --c-bg-soft: #f6f6f7;
  --c-text: #213547;
  --c-text-2: #666;
  --c-brand: #3451b2;
  --c-divider: #e2e2e3;
  --c-ac: #22863a;
  --c-fail: #cb2431;
  --c-warn: #b08800;
  --c-gray: #6a737d;
}

@media (prefers-color-scheme: dark) {
  :root {
    --c-bg: #1a1a1a;
    --c-bg-soft: #252529;
    --c-text: #ddd;
    --c-text-2: #aaa;
    --c-brand: #6a9fff;
    --c-divider: #3a3a3c;
    --c-ac: #56d364;
    --c-fail: #f85149;
    --c-warn: #d29922;
    --c-gray: #8b949e;
  }
}

* { margin: 0; padding: 0; box-sizing: border-box; }

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
  color: var(--c-text);
  background: var(--c-bg);
  line-height: 1.7;
}

a { color: var(--c-brand); text-decoration: none; }
a:hover { text-decoration: underline; }

.top-nav {
  display: flex;
  align-items: center;
  gap: 1.5rem;
  padding: 0.75rem 1.5rem;
  border-bottom: 1px solid var(--c-divider);
  background: var(--c-bg);
  position: sticky;
  top: 0;
  z-index: 10;
}

.site-title { font-weight: 700; font-size: 1.1rem; white-space: nowrap; }

#search { flex: 1; max-width: 400px; position: relative; }
#search .pagefind-ui__form { position: relative; }
#search .pagefind-ui__drawer {
  position: absolute; top: 100%; left: 0; z-index: 100;
  background: var(--c-bg); border: 1px solid var(--c-divider); border-radius: 6px;
  box-shadow: 0 4px 16px rgba(0,0,0,0.15); max-height: 70vh; overflow-y: auto;
  margin-top: 0.25rem; width: max(100%, 500px);
}
#search .pagefind-ui__message, #search .pagefind-ui__results { padding: 0.5rem; font-size: 0.85rem; }
#search .pagefind-ui__result { padding: 0.5rem 0; border-bottom: 1px solid var(--c-divider); }
#search .pagefind-ui__result-title { font-size: 0.9rem; }
#search .pagefind-ui__result-excerpt { font-size: 0.8rem; color: var(--c-text-2); }

.layout {
  display: flex;
  max-width: 1400px;
  margin: 0 auto;
}

.sidebar {
  width: 280px;
  flex-shrink: 0;
  padding: 1rem;
  border-right: 1px solid var(--c-divider);
  height: calc(100vh - 49px);
  overflow-y: auto;
  position: sticky;
  top: 49px;
  font-size: 0.9rem;
}

.sidebar ul { list-style: none; padding-left: 1rem; }
.sidebar > ul { padding-left: 0; }
.sidebar li { margin: 0.15rem 0; }
.sidebar a { display: block; padding: 0.15rem 0; }
.sidebar summary { cursor: pointer; font-weight: 600; padding: 0.15rem 0; }
.sidebar-dir { color: var(--c-text-2); font-size: 0.75em; font-weight: 400; margin-left: 0.5em; }

.content {
  flex: 1;
  min-width: 0;
  padding: 2rem;
}

.content h1 { font-size: 1.8rem; margin-bottom: 1rem; border-bottom: 1px solid var(--c-divider); padding-bottom: 0.5rem; }
.content h2 { font-size: 1.4rem; margin: 2rem 0 0.75rem; }
.content p { margin: 0.5rem 0; }
.content ul, .content ol { margin: 0.5rem 0; padding-left: 1.5rem; }
.content pre { margin: 1rem 0; border-radius: 6px; overflow-x: auto; }
.content code { font-family: 'SF Mono', Monaco, Consolas, monospace; font-size: 0.9em; }
.content img { max-width: 100%; }
.content table { border-collapse: collapse; margin: 0.5rem 0; font-size: 0.9rem; }
.content th, .content td { border: 1px solid var(--c-divider); padding: 0.4rem 0.75rem; }
.content th { background: var(--c-bg-soft); }

.table-wrapper { overflow-x: auto; margin: 0.5rem 0; }

.verify-matrix { font-size: 0.75rem; }
.verify-matrix td, .verify-matrix th { white-space: nowrap; padding: 0.25rem 0.5rem; }

.status-ac { color: var(--c-ac); }
.status-fail { color: var(--c-fail); }
.status-warn { color: var(--c-warn); }
.status-gray { color: var(--c-gray); }

.dot { font-size: 0.6em; vertical-align: middle; }
.dot-ac { color: var(--c-ac); }
.dot-fail { color: var(--c-fail); }
.dot-warn { color: var(--c-warn); }
.dot-gray { color: var(--c-gray); }


.hero { text-align: center; padding: 4rem 1rem; }
.hero h1 { font-size: 2.5rem; border: none; }
.hero p { color: var(--c-text-2); font-size: 1.2rem; margin: 1rem 0 2rem; }
.hero-actions { display: flex; gap: 1rem; justify-content: center; }
.btn { padding: 0.6rem 1.5rem; border-radius: 6px; font-weight: 600; border: 1px solid var(--c-divider); }
.btn-primary { background: var(--c-brand); color: #fff; border-color: var(--c-brand); }

details { margin: 0.5rem 0; }
summary { cursor: pointer; }

.dep-list { list-style: none; padding-left: 1rem; }
.dep-list li { margin: 0.15rem 0; }
.dep-path { color: var(--c-text-2); font-size: 0.85em; margin-left: 0.5em; }

/* shiki dual theme
   ライト: インラインの color/background-color をそのまま使う
   ダーク: CSS変数 --shiki-dark で上書き */
.shiki-wrapper { position: relative; }
.shiki-wrapper .copy-btn {
  position: absolute; top: 0.5rem; right: 0.5rem;
  padding: 0.25rem 0.5rem; border: 1px solid var(--c-divider); border-radius: 4px;
  background: var(--c-bg); color: var(--c-text-2); cursor: pointer; font-size: 0.75rem;
  opacity: 0; transition: opacity 0.2s;
}
.shiki-wrapper:hover .copy-btn { opacity: 1; }
.bundle-btn {
  display: inline-block; margin-bottom: 0.5rem;
  padding: 0.25rem 0.75rem; border: 1px solid var(--c-divider); border-radius: 4px;
  background: var(--c-bg-soft); color: var(--c-text-2); cursor: pointer; font-size: 0.8rem;
}
.bundle-btn:hover { border-color: var(--c-brand); color: var(--c-brand); }
.shiki { padding: 1rem; border-radius: 6px; overflow-x: auto; font-size: 0.85rem; background: #f6f8fa !important; border: 1px solid var(--c-divider); counter-reset: line; }
.shiki .line { display: inline-block; width: 100%; }
.shiki .line::before {
  content: counter(line); counter-increment: line;
  display: inline-block; width: 2.5em; margin-right: 1em;
  text-align: right; color: var(--c-text-2); opacity: 0.5;
  font-size: 0.8em; user-select: none;
}
@media (prefers-color-scheme: dark) {
  .shiki { background: #24292e !important; border-color: var(--c-divider); }
  .shiki, .shiki span { color: var(--shiki-dark) !important; }
}

@media (max-width: 768px) {
  .sidebar { display: none; }
  .content { padding: 1rem; }
}

/* procon-judge の提出の表 */
.judge-note { color: var(--c-text-2); font-size: 0.85rem; margin: 0.25rem 0 0.5rem; }
.judge-table .judge-stale { color: var(--c-text-2); }
.judge-chip { margin-left: 0.4em; padding: 0 0.4em; border: 1px solid var(--c-divider); border-radius: 8px; font-size: 0.85em; color: var(--c-text-2); }
/* Codeforces 相当の CPU で代わりの経路が通らなかった印。procon-judge の順位表と同じく赤くする。 */
.judge-chip-bad { color: var(--c-fail); border-color: var(--c-fail); }
`;

main().catch(console.error);
